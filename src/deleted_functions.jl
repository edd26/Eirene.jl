
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



##########################################################################################

####	TESTING AND DIAGNOSTICS // FUNCTIONWISE

##########################################################################################

function persistencestats(x)
	L = length(x)
	A = Array{Float64}(undef,8,L)
	for ip = 1:L
		i = x[ip]
		pcloud = rand(20,i)
		res = @timed persistf2vr(pcloud,5,model = "pc",fastop=false,record="all")
  		D = res[1]
   		t = res[2]
		A[1,ip] = length(D["prepairs"][5])
		A[2,ip] = length(D["farfaces"][5])
		A[3,ip] = binom(i-1,4)
		A[4,ip] = length(D["trv"][5])
		A[5,ip] = t
		A[6,ip] = i
		A[7,ip] = length(D["farfaces"][5]) - length(D["prepairs"][5]) / binom(i-1,4)
		A[8,ip] = length(D["trv"][5]) / binom(i-1,4)
	end
	return A
end

function persistencemultistats(x;numtrials = 10)
	trials = Array{Any,1}(numtrials)
	for i = 1:numtrials
		trials[i] = persistencestats(x)
	end
	return trials
end




##########################################################################################

####	BATCH OPERATIONS

##########################################################################################

function eirene_batchcsv(
	inputdirectory,
	outputdirectory;
	maxdim = 1,
	model="dmat",
	entryformat="textfile",
	lowerlim=-Inf,
	upperlim=Inf,
	numrad=Inf,
	fastop=true,
	record="cyclerep",
	pointlabels=[],
	verbose=false)

	filenames = readdir(inputdirectory)

	for i = 2:length(filenames)
		filename = filenames[i]
		filepath = "$(inputdirectory)/$(filename)"
		C = eirene(CSV.read(filepath),
				maxdim 		=maxdim,
				model		=model,
				entryformat	=entryformat,
				lowerlim	=lowerlim,
				upperlim	=upperlim,
				numrad		=numrad,
				fastop		=fastop,
				record		=record,
				pointlabels	=pointlabels,
				verbose		=verbose)
		savepath = "$(outputdirectory)/$(filename).jld"
		JLD.save(savepath,"C",C)
	end
end

function savefacespecial(ct::Array{Int64,1},kk::Int64,colsum::Array{Int64,1},farfilt::Int64,oldclaw::Array{Int64,1},rt::Array{Int64,1},zt::Array{Int64,1})
	keep = true
	for l = ct[kk]:colsum[kk]
		if  zt[l]>= farfilt && oldclaw[rt[l]]>=farfilt
			keep = false
			println(["rt[l]" rt[l]])
			break
		end
	end
	if keep
		println("kept")
	end
	return keep
end


function sparsifydesparsifytest(m,n)
# 	added 12/27/2017
	for p = 1:m
		A = rand(n,n).<0.1
		A = convert(Array{Int64},A)
		rv,cp = full2ss(A)
		B = ss2full(rv,cp,n)
		if A != B
			print("error on iteration $(p)")
			return A,B
		end
	end
	print("test successful")
end

function ss2full(rowval,colptr,m)
# 	renamed from 'showfull' on 12/27/2017; ss stands for 'sparse support'
	n = length(colptr)-1
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function ss2full(rowval,colptr,m,n)
# 	renamed from 'showfull' on 12/27/2017; ss stands for 'sparse support'
	M = zeros(Int8,m,n)
	for j = 1:n
		M[rowval[cran(colptr,j)],j]=1
	end
	return M
end

function binom_float(x,y)
	k = 1;
	a = x
	b = 1
	for i = 1:y
		k = k*a/b
		a-=1
		b+=1
	end
	return k
end

function undercat(X)
	l = length(X);
	if l <= 1
		return X
	else
		m = l+l-1;
		Y = Array{Any}(undef,m)
		Y[2:2:m] = "_"
		Y[1:2:m] = X
		return string(Y...)
	end
end

# stands for number of simplicies of cardinality less than or equal to k
# A is an array of arrays
# k is an integer
# this function is not used at the time of this writing (jan 14, 2018)
function numsimcardlek(A,k)
	c = 0;
	for p = 1:k
		c 	+= length(A[p])
	end
	return c
end


function csvimport2linends(M)
	m,n = size(M)
	endpoints = zeros(Int64,m)
	for p = 1:m
		for q = 1:n
			if M[p,q] == ""
				endpoints[p] = q-1
				break
			elseif q == n
				endpoints[p] = n
			end
		end
	end
	return endpoints
end

function printsize(var,varname)
	println(string("size(",varname,") = ",size(var)))
end

function cyclevertices(
	D::Dict;
	dim = 1,
	cycle = 1)

	sd 			= dim+2
	rep 		= getcycle(D,sd,cycle)
	vertices 	= incidentverts(D::Dict,sd-1,rep)
	vertices 	= D["nvl2ovl"][vertices]
	return vertices
end


function getrepsize(D::Dict,classnumber;dim=1)
	sd = dim+2
	if !haskey(D,"cyclerep")
		println("This object does not contain data about cycle representatives.")
		return
	elseif typeof(classnumber)<: Number
		return length(D["cyclerep"][dim+2][classnumber])
	else
		l = length(classnumber)
		rsize = Array{Int64}(undef,l)
		for i = 1:l
			rsize[i] = length(D["cyclerep"][dim+2][classnumber[i]])
		end
		return rsize
	end
end


function birthtime(C;chain=zeros(Int64,0),dim=1)
	if isempty(chain)
		return -Inf
	else
		translator 	= 	C["ocg2rad"]
		sd 	= 	dim+1
		ocg 		= 	C["grain"][sd] # ocg stands for order-canonical grain
		return 		empteval(maximum,translator[ocg[chain]],-Inf)
	end
end

function deathtime(C;chain=zeros(Int64,0),dim=1)
	if isempty(chain)
		return -Inf
	elseif !isempty(chainboundary(C,chain=chain,dim=dim))
		return Inf
	else
		sd 	= 	dim+1;
		Lrv 		= 	C["Lrv"][sd+1]
		Lcp 		= 	C["Lcp"][sd+1]
		Rrv 		= 	C["Rrv"][sd+1]
		Rcp 		= 	C["Rcp"][sd+1]
		numrows 	= 	boundarycorank(C,dim=dim)
		translator 	= 	zeros(Int64,complexrank(C,dim=dim))
		translator[C["tid"][sd+1]] = 1:numrows
		rv 			= 	findall(!iszero,translator[chain])
		rv 			=  	translator[chain[rv]]
		cp 			= 	[1,length(rv)+1]
		rv,cp 	=  	spmmF2silentLeft(Lrv,Lcp,rv,cp,numrows)
		rv,cp 	=  	spmmF2silentLeft(Rrv,Rcp,rv,cp,numrows)
		#
		# recall that plo = the (ordered) subvector consisting of
		# the first (rank(boundary operator)) elements of tid
		#
		if maximum(rv) > length(C["plo"][sd+1])
			return 		Inf
		else
			ocg2rad = 	C["ocg2rad"]
			grain 	= 	C["grain"][sd+1]
			names 	= 	C["phi"][sd+1]
			return 		empteval(maximum,ocg2rad[grain[names[rv]]],-Inf)
		end
	end
end

function chainboundary(C;chain=zeros(Int64,0),dim=1)
	m 			= 	length(chain)
	if m == 0
		return zeros(Int64,0)
	end
	crv 		= 	convert(Array{Int64,1},1:m);
	ccp 		= 	[1,m+1]
	brv,bcp 	= 	boundarymatrix(C,dim=dim,cols=chain)
	brv,bcp 	= 	spmmF2(brv,bcp,crv,ccp,empteval(maximum,brv,0))
	return 			brv
end

function getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numlows    = length(farfaces[sd-1])
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	numclasses = length(cyclenumber)
	summands = Array{Array{Int64,1},1}(undef,numclasses)
	rep 	 = Array{Int64}(undef,numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
		append!(summands[i],[tid[sd][cyclenumber[i]]])
		summandsupp[summands[i]]=true
	end

	lowgenerators = findall(summandsupp)
	numlowgenerators = length(lowgenerators)
	translator = zeros(Int64,numlows)
	translator[lowgenerators] = 1:length(lowgenerators)

	lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

	supp = falses(numlowlows)
	m = size(lowfacemat,1)
	plow2phigtranslator = Array{Int64}(undef,numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	for i = 1:numclasses

		supp[:].= false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end

		brv = findall(supp[tid[sd-1]])
		bcp = [1,length(brv)+1]
		brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
		brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

		rep[i] = length(brv)+length(summands[i])
	end
	return rep
end

function getcyclesize(D::Dict,cyclenumber;dim = 1)
	sd = dim+2
	if !haskey(D,"Lirv")
		if !haskey(D,"cyclerep")
			println("This object does not store a complete cycle basis.")
		else
			println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
		end
		return
	end
	farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]

	rrv = getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	return rrv
end


function nnzbars_test()
	for p = 1:20
		C 	= 	eirene(rand(20,50),model="pc",maxdim=2)
		cr 	= 	C["cyclerep"]
		for q 	=	1:3

		end
	end
end


function integersinsameorderbycolumn(v::Array{Int64,1},maxradue::Int64,colptr)
	# Returns a permutation z on {1,...,length(v)} so that
	# (a) cran(colptr,j) maps to cran(colptr,j) for all j, and
	# (b) crows(colptr,v[z],j) is an array in sorted order
	numcols = length(colptr)-1
	m = length(v)
	x = Array{Int64}(undef,maxradue)
	y = Array{Int64}(undef,maxradue+1)
	z = Array{Int64}(undef,length(v))
	for j = 1:numcols
		x[:] .= 0
		for i = colptr[j]:(colptr[j+1]-1)
			x[v[i]]+=1
		end
		y[1] = colptr[j]
		for i = 1:maxradue
			y[i+1]=y[i]+x[i]
		end
		for i = colptr[j]:(colptr[j+1]-1)
			u = v[i]
			z[i] = y[u]
			y[u]+=1
		end
	end
	return z
end


function intervalcomplementuniquesortedinput(sortedVecOfUniqueIntegers,intervalEndpoint)
	v = sortedVecOfUniqueIntegers
	n = intervalEndpoint
	L = length(v)
	if L==0
		return 1:n
	elseif L==n
		return Array{Int64}(undef,0)
	else
		boundMarker = 1
		upperBound = v[boundMarker]
		complement = Array{Int64}(undef,n-L)
		marker = 0
		for i = 1:(v[end]-1)
			if i<upperBound
				marker+=1
				complement[marker] = i
			else
				boundMarker+=1
				upperBound = v[boundMarker]
			end
		end
		complement[(marker+1):end]=(v[end]+1):n
	end
	return complement
end


function getstarweights(symmat)
	m = size(symmat,1)
	w = zeros(Int64,m)
	getstartweights_subr2(symmat::Array{Int64,2},w::Array{Int64,1},m::Int64)
	return w
end

function getstartweights_subr2(symmat::Array{Int64,2},w::Array{Int64,1},m::Int64)
	s = copy(symmat)
	for i = 1:m
		s[i,i]=0
	end
	l = Array{Int64}(undef,m)
	lDown = Array{Int64}(undef,m)
	val 		= Array{Array{Int64,1}}(undef,m)
	supp 		= Array{Array{Int64,1}}(undef,m)
	suppDown 	= Array{Array{Int64,1}}(undef,m)
	for i = 1:m
		supp[i] = findall(!iszero,s[:,i])
		suppDown[i] = i .+findall(!iszero,s[(i+1):end,i])
		val[i] = s[supp[i],i]
		l[i] = length(supp[i])
		lDown[i] = length(suppDown[i])
	end

	for i = 1:m
		Si = supp[i]
		Vi = val[i]
		for jp = 1:l[i]
			j = Si[jp]
			dij = Vi[jp]
			Sj = suppDown[j]
			Vj = val[j]
			for kp = 1:lDown[j]
				k = Sj[kp]
				if k == i
					continue
				end
				dkj = Vj[kp]
				dki = s[i,k]
				if dki >= dkj && dij >= dkj
					w[i]+=1
				end
			end
		end
	end
	return w
end

function ordercanonicalform(
	S::Array{Tv};
	minrad=-Inf,
	maxrad=Inf,
	numrad=Inf,
	fastop::Bool=true,
	verbose::Bool=false) where Tv

    symmat_float = convert(Array{Float64},copy(S))
    symmat = copy(S)
	m = size(symmat,1);
	convert(Tv,minrad)
	convert(Tv,maxrad)

	effectivemin = -maxrad
	effectivemax = -minrad
	symmat = -symmat

	for i = 1:m
		symmat[i,i]=Inf
	end
	if fastop
		maxmin = -Inf
		for j=1:m
			holdmin = minimum(symmat[:,j])
			if holdmin > maxmin
				maxmin = holdmin
			end
		end
		if maxmin > effectivemin
			effectivemin = maxmin
		end
	end
	if numrad == 1
		for i = 1:m
			for j = (i+1):m
				sij = symmat[i,j]
				if sij>effectivemin
					symmat[i,j]=1
					symmat[j,i]=1
				else
					symmat[i,j]= 0
					symmat[j,i]= 0
				end
			end
		end
		for i = 1:m
			symmat[i,i]= 0
		end
		ocg2rad = [1]
		return round.(Int64,symmat),ocg2rad
	end
	numfilt = binom(m,2)
	for i = 1:m
		for j = (i+1):m
			sij = symmat[i,j]
			if sij >= effectivemax
				symmat[i,j] = Inf
				symmat[j,i] = Inf
			elseif sij < effectivemin
				# note that loose inequality
				# here could cause some bars
				# that disappear at the last
				# grain to appear to
				# live forever
				symmat[i,j] = -Inf
				symmat[j,i] = -Inf
			end
		end
	end
	for i = 1:m
		symmat[i,i] = -Inf
	end
	ocg2rad = zeros(Float64,binom(m,2))
	p = sortperm(symmat[:],alg=MergeSort)
	if verbose
		print("done sorting")
	end
	ordervalue = 0
	floatingvalue = symmat[p[1]]
	for i = 1:m^2
		ii = p[i]
		if symmat[ii] == floatingvalue
			symmat[ii] = ordervalue
		else
			ordervalue+=1
			floatingvalue = symmat[ii]
			ocg2rad[ordervalue]=symmat_float[ii]
			symmat[ii]=ordervalue
		end
	end
	deleteat!(ocg2rad,(ordervalue+1):binom(m,2))
	return round.(Int64,symmat),ocg2rad
end


# Under development as of 12/30/2017
#
# function 	graduate(
# 			A;
# 			minval		=	minimum(A),
# 			maxval		= 	Inf,
# 			numval		= 	Inf,
# 			stepsize	= 	[],
# 			privatemax	= 	Inf)
#
# 	if numval == 1
# 		return ones(Int64,size(A)...),minimum(A)
# 	end
#
# 	#	Note that maxval and numval cannot both be specified by the user
# 	if stepsize 	= 	[]
# 			alpha	=	(maxval-minval)/(numval-1)
# 		end
# 	else
# 		alpha 		= 	stepsize
# 		if  (maxval 	== 	Inf) 	||	(numval == [])
# 			numval	= 	Inf
# 		else
# 			maxval 	= 	minval+(numval-1)*alpha
# 		end
# 	end
#
# 	p 						= sortperm(vec(A))
#
# 	# Compute the ocf
# 	val						= minval
# 	ocg2rad 				= Array{Float64}(undef,size(A,1)*size(A,2)) #the plus 1 covers values taken from the diagonal
# 	ocg2rad[1]				= val
# 	post 					= 1
# 	exceededmax 			= false
# 	ocf						= fill(-1,m,m) #will reverse the order on this after it's been filled
# 	stepcounter = 1
# 	for i = 1:length(p)
# 		if A[p[i]] <= val
# 			ocf[p[i]] = post
# 		else
# 			if numval == Inf
# 				val = A[p[i]]
# 			else
# 				if A[p[i]] == Inf
# 					val = Inf
# 				else
# 					while A[p[i]] > val
# 						stepcounter+=1
# 						if stepcounter == numval
# 							val = maxval # must take this rather cumbersome final step b/c of numerical error, since o/w it can and does happen that the (numval)th grain value fails to equal maxval
# 						else
# 							val+=alpha
# 						end
# 					end
# 				end
# 			end
# 			post+=1
# 			ocf[p[i]] 		= post
# 			ocg2rad[post]	= val
# 		end
# 		if val > privatemax
# 			ocf[p[i:end]] 	= post
# 			exceededmax		= true
# 			break
# 		end
# 	end
# 	if exceededmax
# 		cutoff = post
# 	else
# 		cutoff = post+1
# 	end
# 	deleteat!(ocg2rad,cutoff:length(ocg2rad))
# 	ocg2rad = reverse(ocg2rad,dims=1)
# 	ocf = cutoff - ocf
# end

# NB: Assumes that the input array S has only finite entries.
# NB: The value for keyword <numrad> must be either a positive integer or Inf
function ordercanonicalform_3(
	S::Array{Tv};
	maxrad=Inf,
	minrad=-Inf,
	numrad=Inf,
	vscale="default",
	fastop::Bool=true,
	verbose::Bool=false) where Tv

    if size(S,1) == 0
    	ocf 		= 	zeros(Int64,0,0)
    	ocg2rad		= 	zeros(Float64,0)
    	return ocf,ocg2rad
    end

    # Format input
    symmat		= convert(Array{Float64},copy(S))
	m 			= size(symmat,1)

	if vscale == "default"
		for i = 1:m
			symmat[i,i] = minimum(symmat[:,i])
		end
	elseif typeof(vscale) <: Array
		vscale = convert(Array{Float64},copy(vscale))
		if length(vscale) != m
			print("Error: keyword <vscale> must take a value of <defualt>, <diagonal>, or <v>, where v is a vector of length equal to the number of vertices in the complex.")
			return
		end
		for i=1:m
			if offdiagmin(symmat,i) < vscale[i]
				print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
				return
			else
				symmat[i,i] = vscale[i]
			end
		end
	elseif vscale 	== 	"diagonal"
		vscale		=	Array{Float64,1}(m)
		for i=1:m
			vscale[i]	=	symmat[i,i]
		end
		# the following is in prnciple unnecessary, but it simplifies rounding
		for i=1:m
			if offdiagmin(symmat,i) < vscale[i]
				print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
				return
			end
		end
	end

	# Deterime the public maxrad
	if maxrad == Inf
		publicmax = maximum(symmat)
	else
		publicmax 	= copy(maxrad)
	end

	# Deterime the public minrad
	publicmin	= minimum(symmat)
	publicmin	= max(publicmin,minrad)

	# It's important that this precede the other cases
	if publicmax < publicmin
		return zeros(Int64,m,m),Array{Float64}(undef,0)
	end

	# This covers all remaining cases where numrad ==1
	if numrad == 1
		privatemax = minimum(maximum(symmat,1))
		ocf = zeros(Int64,m,m)
		if fastop && publicmax >= privatemax
			index = findfirst(maximum(symmat,1),privatemax)
			ocf[index,:]=1
			ocf[:,index]=1
			for i = 1:m
				ocf[i,i]=1
			end
		else
			ocf[symmat.<=publicmax]=1
		end
		return ocf,[publicmax]
	end

	# If necessary, determine step size.  Recall we have already treated every case where numrad == 1.
	if numrad < Inf
		alpha 		= (publicmax-publicmin)/(numrad-1)
	elseif numrad == Inf
		alpha 		= Inf
	end

	# If one stops early, determine when to stop
	if fastop
		privatemax = minimum(maximum(symmat,1))
		privatemax = min(privatemax,publicmax)
		if numrad < Inf
			post = publicmin
			stepcounter = 1
			while post < privatemax
				stepcounter+=1
				if stepcounter == numrad
					post = publicmax # must take this rather cumbersome step on account of numerical error.
				else
					post+=alpha
				end
			end
			privatemax = post
		end
	else
		privatemax = publicmax
	end

	# Extract sortperm
	p 						= sortperm(vec(symmat),alg=MergeSort)

	# Compute the ocf
	val						= publicmin
	ocg2rad 				= Array{Float64}(undef,binom(m,2)+m) #the plus 1 covers values taken from the diagonal
	ocg2rad[1]				= val
	post 					= 1
	exceededmax 			= false
	ocf						= fill(-1,m,m) #will reverse the order on this after it's been filled
	stepcounter = 1
	for i = 1:length(p)
		if symmat[p[i]] <= val
			ocf[p[i]] = post
		else
			if numrad == Inf
				val = symmat[p[i]]
			else
				if symmat[p[i]] == Inf
					val = Inf
				else
					while symmat[p[i]] > val
						stepcounter+=1
						if stepcounter == numrad
							val = publicmax # must take this rather cumbersome final step b/c of numerical error, since o/w it can and does happen that the (numrad)th grain value fails to equal publicmax
						else
							val+=alpha
						end
					end
				end
			end
			post+=1
			ocf[p[i]] 		= post
			ocg2rad[post]	= val
		end
		if val > privatemax
			ocf[p[i:end]] 	= post
			exceededmax		= true
			break
		end
	end
	if exceededmax
		cutoff = post
	else
		cutoff = post+1
	end
	deleteat!(ocg2rad,cutoff:length(ocg2rad))
	ocg2rad = reverse(ocg2rad,dims=1)
	ocf = cutoff - ocf
	return ocf,ocg2rad # additional outputs for diagnostic purposes -> #,privatemax,S,maxrad,publicmax,publicmin,maxrad
end


function truncatearray(N,minrad,maxrad)
	S 					= 	copy(N) #NB it has been verified experimentally that it is VERY important to use the copy function here
	S[S.<minrad] 	   .= 	minrad
	S[S.>maxrad] 	   .= 	Inf
	return 					S
end


#=
- returns finite values for mingrid and maxgrid
- it will always be true that mingrid <= maxgrid
=#
function ceil2grid_overflowparameters(	N;
										minrad 	=	-Inf,
										maxrad	=	Inf
										)
	if (minrad == Inf) || (maxrad == -Inf)
		mingrid 	= 	0
		maxgrid 	= 	0
		return 			minrad,maxrad,mingrid,maxgrid
	end

	S 								= 	Array{Float64}(copy(N)) #NB it has been verified experimentally that it is VERY important to use the copy function here

	if 	minrad 						== 	"minedge"
		minrad 						= 	minimum(offdiagmin(S))
	end

	S[S.<minrad] 				   .= 	minrad
	S[S.>maxrad] 				   .= 	maxrad

	fi 								= 	findall(isfinite,S)
	if 									isempty(fi)
		mingrid 					= 	0
		maxgrid 					= 	0
	else
		# note the case minrad == Inf has already been covered,
		# as has the case where minimum(S[fi]) might be infinite
		if minrad 					==	-Inf
			mingrid 				= 	minimum(S[fi])
		else
			mingrid 				= 	minrad
		end
		# note the case maxrad == -Inf has already been covered
		# as has the case where maximum(S[fi]) might not be infinite
		if maxrad 					== 	Inf
			maxgrid					= 	maximum(S[fi])
		else
			maxgrid					= 	maxrad
		end
	end

	return	minrad, maxrad, mingrid, maxgrid
end

function boundarymatrix(C;dim=1,rows="a",cols="a")
	crr 					= 	complexrank(C,dim=dim-1)
	crc 					= 	complexrank(C,dim=dim)
	if rows == "a"
		rows 				= 	1:crr#complexrank(C,dim=dim-1)
	end
	if cols == "a"
		cols 				= 	1:crc#complexrank(C,dim=dim)
	end
	if empteval(maximum,cols,0) > crc
		println()
		println("error: keyword argument <cols> contains an integer greater than the rank of the complex in dimension <dim>")
		println()
		return
	elseif empteval(maximum,rows,0) > crc
		println()
		print("error: keyword argument <rows> contains an integer greater than the rank of the complex in dimension <dim-1>")
		println()
		return
	end
	if isempty(rows) || isempty(cols)
		rv 					= 	zeros(Int64,0)
		cp 					= 	ones(Int64,length(cols)+1)
		return 					rv,cp
	end
	ncols 					= 	length(cols)
	nrows 					= 	length(rows)
	sd 						= 	dim+1;
	if haskey(C,"farfaces")
		rv 					= 	ff2aflight(C,dim+1,cols)
		rv  				= 	reshape(rv,length(rv))
		cp  				= 	convert(Array{Int64,1},1:sd:(ncols*sd+1))
		cols 				= 	1:ncols
	else
		rv 					= 	C["rv"][sd]
		cp 					= 	C["cp"][sd]
	end
	rv,dummrv,cp,dummycp 	= 	stackedsubmatrices(
								rv,
								cp,
								rows,
								Array{Int64}(undef,0),
								cols,
								max(empteval(maximum,rows,0),empteval(maximum,rv,0))
								)
	return 						rv,cp
end



function boundaryrank(C;dim=1)
	sd 	= 	dim+1;
	if complexrank(C,dim=dim) == 0
		return 0
	else
		return length(C["plo"][sd])
	end
end

function boundarycorank(C;dim=1)
	sd 	= 	dim+1;
	if complexrank(C,dim=dim) == 0
		return 0
	else
		return complexrank(C,dim=dim)-boundaryrank(C,dim=dim)
	end
end

