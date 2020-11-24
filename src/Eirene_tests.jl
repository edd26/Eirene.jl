# created 2018-04-15
function unittest()

	numits 	= 	5
	maxdim 	= 	2
	x 		= 	Array{Any}(undef,24)

	x[1] 	= 	eirenevrVperseusvr() 					# correct answer: empty
	x[2] 	= 	eirenevrVeirenepc(numits,maxdim) 		# correct answer: empty
	x[3] 	= 	eirenevrVeirenecomplex(numits,maxdim)	# correct answer: empty
	x[4] 	= 	eirenecomplexVhandcalc() 				# correct answer: empty
	x[5] 	= 	eirenevrVeirenesuspendedcomplex(numits) # correct answer: empty
	x[6] 	= 	checkparameters()						# correct answer: empty
	x[7] 	= 	checkcomplexformats()					# correct answer: empty
	x[8] 	= 	checksuspension(numits)					# correct answer: empty
	x[9] 	= 	checkoffdiagmin(numits)					# correct answer: empty
	x[10]	= 	checkoffdiagmean(numits)				# correct answer: empty
	x[11] 	= 	checkceil2grid(numits)					# correct answer: empty
	x[12]	= 	checkminmaxceil(numits) 				# correct answer: empty
	x[13] 	= 	checksegmentarray(numits)				# correct answer: empty
	x[14]	= 	checksegVdesegcomplex(numits)			# correct answer: empty
	x[15]	= 	checkeulervector2dimensionpattern()		# correct answer: empty
	x[16]	= 	checkdimensionvalues2dimensionpattern() # correct answer: empty
	x[17]	= 	checkbuildcomplex3_diagentries(numits) 	# correct answer: empty
	x[18] 	= 	checktrueordercanonicalform(numits) 	# correct answer: empty
	x[19]	= 	checkloadfile()							# correct answer: empty
	x[20]   =   wd_test_1()								# correct answer: empty
	x[21]   =   wd_test_2()								# correct answer: empty
	x[22]   =   wd_test_3()								# correct answer: empty
	x[23] 	= 	wd_test_4()								# correct answer: empty
	x[24]	=	wd_test_5()								# correct answer: empty

	for p 	= 	1:length(x)
		if !isempty(x[p])
			println(p)
			return x
		end
	end
	return []
end

#1 
function eirenevrVperseusvr()
	passtest 			= 	true

	datapath 			=   Eirene.testfp("prsjd")
	E 					= 	JLD.load(datapath)
	E 					= 	E["E"]

	for p 				= 	1:length(E)
		vrmatr 			= 	E[p]["vrmat"]
		maxdim 			= 	E[p]["maxdim"]-1	# for perseus, maxdim refers to simplices, not homlogy
		prdict 			= 	E[p]["perseusdict"] 						# persesus dictionary
		endict 			= 	Eirene.eirene(vrmatr,model="vr",maxdim=maxdim)		# eirene dictionary

		i,j 			= 	Eirene.firstbcdiff([endict prdict],maxdim=maxdim)
		if 	i 			!= 	0
			return 		endict, prdict, j
		end
	end
	return 				[]
end

#2
function eirenevrVeirenepc(numits,maxdim)
	for p 			= 	1:numits
		numpts 		= 	rand(50:100,1)
		numpts 		= 	numpts[1]
		ambdim 		= 	rand(1:100,1)
		ambdim 		= 	ambdim[1]
		pc 			= 	rand(ambdim,numpts)
		d 			= 	colwised(pc)
		nodrad		= 	offdiagmin(d)./2

		Cvr 		= 	eirene( d,model="vr",maxdim=maxdim)
		Cpc 		= 	eirene(pc,model="pc",maxdim=maxdim,nodrad=nodrad)

		i,j 		= 	firstbcdiff([Cpc Cvr],maxdim=maxdim)
		if 	i 		!= 	0
			return  Cpc, Cvr, j
			break
		end
	end
	return 				[]
end

#3
function eirenevrVeirenecomplex(numits,maxdim,verbose=false)
	for p 			= 	1:numits
		numpts 		= 	rand(50:60,1)
		numpts 		= 	numpts[1]
		d 			= 	vertexlifemat(numpts,model="rand")
		Cvr 		= 	eirene(d,maxdim=maxdim,model="vr")

		if verbose
			println("eirenevrVeirenecomplex - waypoint 1")
		end

		rv,cp 		= 	boundarymatrices(Cvr)
		fv 			= 	ocff2of(Cvr["grain"],Cvr["ocg2rad"])

		if verbose
			println("eirenevrVeirenecomplex - waypoint 2")
		end

		Ccx		 	= 	eirene(rv=rv,cp=cp,fv=fv,model = "complex",maxdim=maxdim)

		if verbose
			println("eirenevrVeirenecomplex - waypoint 3")
		end

		i,j 		= 	firstbcdiff([Cvr Ccx],maxdim=maxdim)

		if verbose
			println("eirenevrVeirenecomplex - waypoint 4")
		end

		if 	i 			!= 	0
			return 	 	Cvr, Ccx, i
		end
	end
	return 				[]
end

#4 
function eirenecomplexVhandcalc()
	passtest 				= 	true
	numits 					= 	2
	maxdim 					= 	2
	maxrad 					= 	100
	numradset 				= 	Array{Any}(undef,2)  # must take extra care to ensure that
	numradset[1]			= 	10
	numradset[2]			= 	Inf
	for space 				= 	["sphere", "empty","skrabatriangle"]
		for numrad 			= 	numradset
			if 		space 	== 	"sphere"
				pathkey 	= 	"hsphr"
				entryformat = 	"dp"
			elseif 	space 	== 	"empty"
				pathkey 	= 	"hempt"
				entryformat = 	"dp"
			elseif space 	== 	"skrabatriangle"
				pathkey 	= 	"hstri"
				entryformat = 	"sp"
			end

			datapath 		= 	testfp(pathkey)

			C 				= 	eirene(	datapath,
										model 		= 	"complex",
										entryformat = 	entryformat,
										maxdim 		= 	maxdim)

			D				= 	handcalcsolution() #  load(solnpath)
			solkey	 		= 	solutionkey(
								model 		= 	"complex",
								maxrad 		= 	maxrad,
								numrad		= 	numrad,
								space		= 	space,
								problemset 	= 	"hand")
			D				= 	D[solkey]

			i,j 			= 	firstbcdiff([C D],maxdim=maxdim)
			if 	i 		   != 	0
				return 			C, D, j
			end
		end
	end
	return					[]
end

#5
function eirenevrVeirenesuspendedcomplex(numits)
	maxdim 			= 	2
	for p 			= 	1:numits
		numpts 		= 	rand(50:60,1)
		numpts 		= 	numpts[1]
		d 			= 	vertexlifemat(numpts,model="rand")
		Cvr	 		= 	eirene(d,maxdim=maxdim,model="vr")

		rv,cp,fv 	= 	eirened2complex(Cvr)
		rvs,cps,fvs = 	suspend(rv,cp,fv,degree=2)
		Csx 		= 	eirene(rv=rvs,cp=cps,fv=fvs,maxdim=maxdim+2,model="complex")

		i,j 		= 	firstbcdiff([Cvr Csx],maxdim=maxdim,offset=2)
		if 	i 			!= 	0
			return 	 	Cvr, Csx, i
		end
	end
	return 				[]
end


#6
function 	checkparameters()
	numpts 		= 	rand(20:60,1)
	numpts 		= 	numpts[1]
	ambdim 		= 	rand(1:100,1)
	ambdim 		= 	ambdim[1]
	pc 			= 	rand(ambdim,numpts)
	d 			= 	colwised(pc)
	nodrad		= 	offdiagmin(d)./2

	C0 			= 	eirene(d,model="vr",maxdim=2)

	for maxdim 		= [0 1]
		for minrad 		= [-Inf 0 1 Inf]
			for maxrad 		= [-Inf 0 1 Inf]
				for fastop 		= [true,false]
					for vscale 		= [[]]
						for record 		=  ["all" "cyclerep" "none"]
							for pointlabels	= [[]]
								#######  add numrad iterable here ######
								# for numrad =

								Cvr  	= 	eirene(
											d;
											model		= "vr",
											maxdim 		= maxdim,
											minrad		= minrad,
											maxrad		= maxrad,
											fastop		= fastop,
											vscale		= vscale,
											record		= record,
											pointlabels	= pointlabels)

								Cpc  	= 	eirene(
											pc;
											model		= "pc",
											maxdim 		= maxdim,
											minrad		= minrad,
											maxrad		= maxrad,
											fastop		= fastop,
											vscale		= vscale,
											record		= record,
											pointlabels	= pointlabels)

								rv,cp 	= 	boundarymatrices(C0)
								fv 		= 	ocff2of(
											C0["grain"],
											C0["ocg2rad"])

								Ccx 	= 	eirene(
											rv 			= 	rv,
											cp			= 	cp,
											fv			= 	fv,
											model 		= 	"complex",
											record 		= 	record)

								if 	minrad 			   == 	-Inf && maxrad 	==	Inf
									X 					= 	[C0 Cvr Cpc Ccx]
								else
									X 					= 	[Cvr Cpc Ccx]
								end
								i,j 					= 	firstbcdiff(X,maxdim=maxdim)
								if i 				   !=	0
									return 					X[1], X[i], j
								end

								if record 				== 	"all"
									for 	dim 		= 	0:maxdim
										if 	!generatorbdc(Cpc,dim=dim)
											return Cpc
										end
										if 	!generatorbdc(Cvr,dim=dim)
											return Cvr
										end
										if 	!generatorbdc(Ccx,dim=dim)
											return Ccx
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
	return	[]
end

#7
function checkcomplexformats()
	pc 						= 	rand(20,60)
	maxdim 					= 	2
	C 						= 	eirene(pc,model="pc",maxdim=maxdim,record="all")

	# format 1
	rvsg,cpsg 				= 	boundarymatrices(C) # sg stands for "segemnted"
	fvsg 					= 	ocff2of(C["grain"],C["ocg2rad"]) # sg stands for "segemnted"

	# format 2
	rv,cp,fv,dp 			=	segmentedfilteredcomplex2unsegmentedfilteredcomplex(rvsg,cpsg,fvsg)

	# format 3
	dv 						=	dimensionpattern2dimensionvalues(dp)

	# format 4
	ev 						=	diff(dp)

	Csg 					= 	eirene(rv=rvsg,cp=cpsg,fv=fvsg,model="complex",record="all")
	Cdp 					= 	eirene(rv=rv,cp=cp,fv=fv,dp=dp,model="complex",record="all")
	Cdv						= 	eirene(rv=rv,cp=cp,fv=fv,dv=dv,model="complex",record="all")
	Cev 					= 	eirene(rv=rv,cp=cp,fv=fv,ev=ev,model="complex",record="all")

	X 						= 	[Csg,Cdp,Cdv,Cev]
	i,j 					= 	firstbcdiff(X,maxdim=maxdim)
	if i 				   	!=	0
		return 					X[1], X[i], j
	end
	passedtest 				= 	true
	for 	dim 		= 	0:maxdim
		if 	!generatorbdc(C,dim=dim)
			passedtest 	= 	false
			break
		end
		if 	!generatorbdc(Csg,dim=dim)
			passedtest 	= 	false
			break
		end
		if 	!generatorbdc(Cdp,dim=dim)
			passedtest 	= 	false
			break
		end
		if 	!generatorbdc(Cdv,dim=dim)
			passedtest 	= 	false
			break
		end
		if 	!generatorbdc(Cev,dim=dim)
			passedtest 	= 	false
			break
		end
	end
	if passedtest
		return zeros(Int64,0)
	else
		return false
	end
end

#8
function checksuspension(numits)
	for p 				= 	1:numits
		for degree 		= 	[0,1,5]
			x 				= 	rand(20,50)
			delrange 		= 	1:degree

			C 				= 	eirene(x,model="pc",maxdim=2)
			rv,cp,fv 		= 	eirened2complex(C)
			C2 				= 	eirene(rv=rv,cp=cp,fv=fv)
			rv2,cp2,fv2 	= 	eirened2complex(C2)
			rv3,cp3,fv3 	= 	suspend(rv2,cp2,fv2,degree=degree)

			check1 			= 	rv3[delrange] == fill(zeros(Int64,0),degree)
			check2 			= 	cp3[delrange] == fill([1],degree)
			check3 			= 	fv3[delrange] == fill(zeros(Int64,0),degree)

			check4 			= 	rv3[degree+1:end] 	== 	rv2
			check5 			= 	cp3[degree+1:end] 	== 	cp2
			check6 			= 	fv3[degree+1:end] 	== 	fv2

			if !all([check1,check2,check3,check4,check5,check6])
				return x,C,rv,cp,fv,C2,rv2,cp2,fv2,rv3,cp3,fv3
			end
		end
	end
	return zeros(Int64,0)
end

#9
function checkoffdiagmin(numits)
	for p 					= 	1:numits
		S 					= 	rand(50,50)
		for q 				= 	1:50
			checkval 		= 	minimum(deleteat!(S[:,q],q))
			if checkval 	!= 	offdiagmin(S,q)
				return S
			end
		end
	end
	return zeros(Int64,0)
end

#10
function checkoffdiagmean(numits)
	for p = 1:numits
		numpts 	= 	rand(50:100,1)
		numpts 	= 	numpts[1]
		S 		= 	rand(numpts,numpts)
		T 		= 	copy(S)
		mu 		= 	offdiagmean(T,defaultvalue = 0)
		u 		= 	zeros(numpts)
		for i 	= 	1:numpts
			ran =   setdiff(1:numpts,i)
			u[i]= 	mean(S[ran,i])
		end
		if 			u[:] != mu[:]
			println()
			println("error: please check <offdiagmean>")
			return S,T,mu,u
		end
	end
	return []
end

#11
function checkceil2grid(numits)
	for 	p 	=	1:numits
		A 		= 	rand(70)
		ran 	= 	sort(rand(20))
		ran 	= 	ran*(1/maximum(ran)) # this guarantees that maximum(ran) > maximum(A)
		crct 	= 	crosscheckceil2grid(A,ran)
		if !crct
			println()
			println("error: please check checkceil2grid")
		end
	end
	return []
end

#12
function checkminmaxceil(numits)
	for p 			= 	1:numits
		numpts 		= 	rand(50:100,1)
		numpts 		= 	numpts[1]
		x 			= 	rand(numpts,numpts)
		x 			= 	x+x';
		numradA		= 	Array{Any}(rand(10:50,2))
		maxradA		= 	Array{Any}(rand(2))
		minradA		= 	Array{Any}(rand(2))
		append!(maxradA,[Inf])
		append!(minradA,[-Inf])
		append!(minradA,["minedge"])
		append!(numradA,[1,Inf])
		for	maxrad in maxradA
			for minrad in minradA
				for numrad in numradA
					A 	= 	minmaxceil(x,minrad=minrad,maxrad=maxrad,numrad=numrad)
					B   = 	ceilvr(x,minrad=minrad,maxrad=maxrad,numrad=numrad)
					if  	A != 	B
						println()
						println("error: please check minmaxceil and ceilvr")
						return x,minrad,maxrad,numrad
					end
				end
			end
		end
	end
	return []
end

#13
function segmentarray(vr::Tv,vp) where Tv
	m 	= 	length(vp)-1
	u 	= 	Array{Tv}(undef,m)
	for 	p 	=	1:m
		u[p] 	= 	crows(vp,vr,p)
	end
	return u
end

function checksegmentarray(numits)
	for p 			= 	1:numits
		cp 			= 	rand(1:100,50)
		dp 			= 	eulervector2dimensionpattern(cp)
		rv 			= 	rand(dp[end]-1)
		A 			= 	segmentarray(rv,dp)
		for q 		= 	1:50
			if 	A[q]	!=	crows(dp,rv,q)
				return rv,cp
			end
		end
	end
	return []
end

#14
function checksegVdesegcomplex(numits)
	for p 			= 	1:numits
		x 			= 	rand(30,30)
		x 			= 	x+x'
		C 			= 	eirene(x,model="vr",maxdim=2)
		rv,cp 		= 	boundarymatrices(C)
		fv 			= 	ocff2of(C["grain"],C["ocg2rad"])
		rvold 		= 	copy(rv)


		rv1,cp1,fv1,dp1 	= 	segmentedfilteredcomplex2unsegmentedfilteredcomplex(rv,cp,fv)
		rv2,cp2,fv2,dp2 	= 	unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv1,cp1,fv1,dp1)

		check0 	= 	dp1 == dp2
		check1 	= 	rv2 == rv
		check2 	= 	cp2 == cp
		check3 	= 	fv2 == fv
		check4 	= 	true
		for q 	= 	1:length(dp1)-1
			if 	length(fv[q]) != length(cran(dp1,q))
				check4 	= 	false
				break
			end
		end

		if !all([check0, check1, check2,check3,check4])
			println([check0, check1,check2,check3,check4])
			return x,rv1,cp1,fv1,dp1,rv2,cp2,fv2,rv,cp,fv,rvold
		end
	end
	return []
end

#15
function checkeulervector2dimensionpattern()
	for p = 1:1000
		l 	= 	rand(100:500,1)
		l 	= 	l[1]
		ev 	= 	rand(1:1000,l)
		dp 	= 	eulervector2dimensionpattern(ev)
		if (dp[1]!=1) || (diff(dp) != ev)
			return ev,dp
		end
	end
	return []
end


#16
function checkdimensionvalues2dimensionpattern()
	numits 			= 	10
	for p 			= 	1:numits
		dv 			= 	zeros(Int64,0)
		for q 		= 	0:100
			x 		= 	rand(0:2)
			x 		= 	x[1]
			append!(dv,q*ones(Int64,x))
		end
		dp1			= 	dimensionvalues2dimensionpattern(dv)
		dv2 		= 	dimensionpattern2dimensionvalues(dp1)
		dp3 		= 	dimensionvalues2dimensionpattern(dv2)
		dv4 		= 	dimensionpattern2dimensionvalues(dp3)

		check1 	= 	(dp1 == dp3)
		check2 	= 	(dv2 == dv4) && (dv2 == dv)

		if !all([check1 check2])
			println("error: please check <dimensionpattern2dimensionvalues>")
			return dv
		end
	end
	return []
end



#17

#=
This function was used to check that
(1) in dimensions higher than 0,
the farfaces, firstv, filtration, prepairs (plus the vertex permutation)
generated by buildcomplex3 are invariant under changes to the diagonal.
(2) the filtration in dimension 0 is an approrpriately permuted copy of the
diagonal of the input matrix
(3) the value assigned to key "symmat" is an approrpriately permuted copy of
the input matrix itself
=#
function checkbuildcomplex3_diagentries(numits)
	for p 				= 	1:numits
		for maxsd 	= 	1:4
			numpts 			= 	rand(10:50,1)
			numpts 			= 	numpts[1]
			d 				= 	rand(1:10^6,numpts,numpts)
			d 				= 	d+d'
			randv 			= 	rand(1:numpts,1)
			randv 			= 	randv[1]
			d[randv,randv] 	= 	maximum(d)+1
			t 				= 	copy(d)
			for 	q 		= 	1:numpts
				t[q,q] 		= 	0
			end

			diagd 			= 	diag(d)
			diagt 			= 	diag(t)

			if any(diag(t).!=0)
				println("please check that the diagonal entries of s are zeros")
			end
			if count(!iszero,t)!= numpts^2-numpts
				println("please check that the off-diagonal entries of s are nonzero")
			end

			Dd 				= 	buildcomplex3(d,maxsd)
			Dt 				= 	buildcomplex3(t,maxsd)

			Ddc 			= 	copy(Dd)
			Dtc 			= 	copy(Dt)
			Ddc["grain"] 	= 	Ddc["grain"][2:end]
			Dtc["grain"] 	= 	Dtc["grain"][2:end]
			Ddc["symmat"] 	= 	[]
			Dtc["symmat"]	= 	[]

			check1 			= 	Ddc == Dtc
			check2 			= 	Dd["grain"][1] 	== 	diagd[Dd["nvl2ovl"]]
			check3 			= 	Dt["grain"][1]	== 	diagt[Dt["nvl2ovl"]]
			check4 			= 	Dd["symmat"] 	== 	d[Dd["nvl2ovl"],Dd["nvl2ovl"]]
			check5 			= 	Dt["symmat"] 	== 	t[Dt["nvl2ovl"],Dt["nvl2ovl"]]

			if !(check1 && check2 && check3 && check4 && check5)
				return Dd,Dt
			end
		end
	end
	return []
end


#18
function checktrueordercanonicalform(numits)
	for p 	= 	1:numits
		m 	= 	rand(50:300,1)
		m 	= 	m[1]
		x 	= 	rand(m,m)
		if 		isodd(p)
			x 	= 	x+x';
		end
		ocf,val 	= 	trueordercanonicalform(x,factor=true)
		check1 		= 	x == val[ocf]
		check2 		= 	val == sort(unique(x))

		k 			= 	rand(1:100,1)
		k 			= 	k[1]
		ocau,valu	=   trueordercanonicalform(x,firstval=1,factor=true,rev=false)
		ocad,vald	=   trueordercanonicalform(x,firstval=1,factor=true,rev=true)
		check3 		= 	ocau == (maximum(ocad)+1).-ocad
		check4 		= 	x == valu[ocau]
		check5 		= 	x == vald[ocad]
		check6 		= 	length(valu) == maximum(ocau)
		check7 		= 	length(vald) == maximum(ocad)

		if !all([check1, check2, check3, check4, check5, check6, check7])
			println("error: please check trueordercanonicalform")
			println([check1, check2, check3, check4, check5, check6, check7])
			return x
		end
	end
	return []
end


#19
function checkloadfile()

	# FOR REFERENCE
	# "csvdp" => joinpath(@__DIR__,"test/fileload/cell_dp.txt"),
	# "csvdv" => joinpath(@__DIR__,"test/fileload/cell_dv.txt"),
	# "csvev" => joinpath(@__DIR__,"test/fileload/cell_ev.txt"),
	# "csvvr" => joinpath(@__DIR__,"test/fileload/vr.csv"),
	# "csvpc" => joinpath(@__DIR__,"test/fileload/pc.csv"),
	# "txtvr" => joinpath(@__DIR__,"test/fileload/vr.txt"),
	# "txtpc" => joinpath(@__DIR__,"test/fileload/pc.txt"),

	fp_spher			= 	testfp("csvdp")
	fp_csvdp 			= 	testfp("csvdp")
	fp_csvdv 			= 	testfp("csvdv")
	fp_csvev 			= 	testfp("csvev")
	fp_csvsp 			= 	testfp("csvsp")

	fp_csvvr 			= 	testfp("csvvr")
	fp_csvpc 			= 	testfp("csvpc")

	V 					= 	Array{Float64}([0 1 2 3])
	d 					= 	Array{Float64}([0 1 2 3;1 0 1 2; 2 1 0 1; 3 2 1 0])
	h 					= 	x -> delete!(x,"input")

	C 					= 	Array{Any}(undef,5)
	C[1] 				= 	eirene(fp_spher,model="complex",entryformat="dp",maxdim=2,record="all")
	C[2] 				= 	eirene(fp_csvdp,model="complex",entryformat="dp",maxdim=2,record="all")
	C[3] 				= 	eirene(fp_csvdv,model="complex",entryformat="dv",maxdim=2,record="all")
	C[4] 				= 	eirene(fp_csvev,model="complex",entryformat="ev",maxdim=2,record="all")
	C[5] 				= 	eirene(fp_csvsp,model="complex",entryformat="sp",maxdim=2,record="all")

	D 					= 	Array{Any}(undef,4)
	D[1] 				= 	eirene(fp_csvvr,model="vr",maxdim=2,record="all")
	D[2]				= 	eirene(d,model="vr",maxdim=2,record="all")
	D[3] 				= 	eirene(fp_csvpc,model="pc",maxdim=2,record="all")
	D[4] 				= 	eirene(V,model="pc",maxdim=2,record="all")

	for p 				= 	1:5
		if !isassigned(C,p)
			println("C is not completely assigned")
			return C
		end
	end

	for p 				= 	1:4
		if !isassigned(D,p)
			println("D is not completely assigned")
			return D
		end
	end

	checkcx 			= 	pairwiseisequal(copy(C),		under=h)
	checkvr 			= 	pairwiseisequal(copy(D[1:2]),	under=h)
	checkpc 			= 	pairwiseisequal(copy(D[3:4]),	under=h)
	checkbc 			= 	true

	for p 				= 	1:length(C)
		if !isempty(comparebarcodesagainstcomplex(C[p]))
			checkbc 	= 	false
			break
		end
	end

	for p 				= 	1:length(D)
		if !isempty(comparebarcodesagainstcomplex(D[p]))
			checkbc 	= 	false
			break
		end
	end

	check = [checkcx,checkvr,checkpc,checkbc]
	if !all(check)
		return C,D,check
	else
		return []
	end
end


function solutionkey(		;
				model 		= 	"complex",
				maxrad 		= 	Inf,
				numrad		= 	Inf,
				space		= 	"sphere",
				problemset 	= 	"hand")
	if 	problemset 			== 	"hand"
		solkey 		= 	string(
						"mod",model,
						"mar",string(maxrad),
						"nr",string(numrad),
						"it",string(space))
	end
	return solkey
end

function suspend!(rv,cp,fv;degree=1)
	F 	= 	fill(Array{Int64,1}(undef,0),degree)
	v   = 	fill([1],degree)
	prepend!(rv,F)
	prepend!(cp,v)
	prepend!(fv,F)
end

function suspend(rv,cp,fv;degree=1)
	rv 			= 	copy(rv)
	cp 			= 	copy(cp)
	fv 			= 	copy(fv)
					suspend!(rv,cp,fv,degree=degree)
	return 			rv,cp,fv
end


function vertexlifemat(d;model="rand",scale=1/2)
	if 		model 	== 	"pc"
		s 	= 	colwised(d)
	elseif model 	== 	"vr"
		s 	= 	copy(d)
	elseif model 	== 	"rand"
		s 	= 	iudsymmat(d)
	end
	v 		= 	offdiagmin(s)
	if typeof(scale) <: Number
		for	p 	= 	1:size(s,2)
			s[p,p]	= 	v[p]*scale
		end
	elseif 	scale == "rand"
		for p 	= 	1:size(s,2)
			r 		= 	rand(1)
			s[p,p] 	=   r[1]*v[p]
		end
	else
		println()
		println("error: scale must be either a scalar or the string \"rand\"")
		return
 	end
	return 	s
end

function testfp(s)
	D 	= 	Dict(
			"prsip" => joinpath(@__DIR__,"perseus/reservoir/input.txt"),
			"prsop" => joinpath(@__DIR__,"perseus/reservoir/output"),
			"prsjd" => joinpath(@__DIR__,"test/perseus/testdata.jld"), # peresus julian data
			"hantx" => joinpath(@__DIR__,"test/handcalc/testdata.txt"),
			"hanjd" => joinpath(@__DIR__,"test/handcalc/testdata.jld"),
			"hsphr" => joinpath(@__DIR__,"test/handcalc/sphere.csv"),
			"hempt" => joinpath(@__DIR__,"test/handcalc/empty.csv"),
			"hstri" => joinpath(@__DIR__,"test/handcalc/skrabatriangle.csv"),
			"vrmat" => joinpath(@__DIR__,"test/fileload/vrmat.txt"), # vietoris rips
			"csvdp" => joinpath(@__DIR__,"test/fileload/cell_dp.csv"),
			"csvdv" => joinpath(@__DIR__,"test/fileload/cell_dv.csv"),
			"csvev" => joinpath(@__DIR__,"test/fileload/cell_ev.csv"),
			"csvsp" => joinpath(@__DIR__,"test/fileload/cell_sm.csv"),
			"csvvr" => joinpath(@__DIR__,"test/fileload/vr.csv"),
			"csvpc" => joinpath(@__DIR__,"test/fileload/pc.csv"),
			)
	return  D[s]
end


function vrmat(C::Dict)
	if C["input"]["model"] != "vr"
		println()
		println("error: <vrmat> only applies to vietoris-rips complexes")
	end
	nvl2ovl 			= 	C["nvl2ovl"]
	numpts 				= 	length(nvl2ovl)
	ovl2nvl 			= 	Array{Int64}(undef,numpts)
	ovl2nvl[nvl2ovl] 	=   1:numpts
	symmat 				= 	copy(C["symmat"])
	symmat 				= 	symmat[ovl2nvl,ovl2nvl]
	s 					= 	Array{Float64}(undef,symmat)
	for p 				= 	1:length(s)
		if symmat[p]	==   0
			s[p] 		= 	Inf
		else
			s[p] 		= 	C["ocg2rad"][symmat[p]]
		end
	end
	return 				s
end


function iudsymmat(m)
	x	= 	rand(m,m)
	for p 	= 	1:m
		for 	q 	= 	1:p-1
			x[q,p] 	= 	x[p,q]
		end
	end
	return x
end


function 	colwised(x)
	return	Distances.pairwise(Euclidean(),x,dims=2)
end


function comparebarcodesagainstcomplex(C)
	rv,cp 		= 	boundarymatrices(C)
	fv 			= 	ocff2of(C["grain"],C["ocg2rad"])
	maxdim 		= 	C["input"]["maxdim"]
	Ccx		 	= 	eirene(rv=rv,cp=cp,fv=fv,model = "complex",maxdim=maxdim)

	i,j 		= 	firstbcdiff([C Ccx],maxdim=maxdim)
	if 	i 			!= 	0
		return 	 	Cvr, Ccx, i
	else
		return 		[]
	end
end

function firstbcdiff(arrayofdicts;maxdim=1,offset=0) # stands for first barcode difference
	A 					= 	arrayofdicts
	for p 				= 	2:length(A)
		q 				= 	firstbcdiff(A[1],A[p],maxdim=maxdim,offset=offset)
		if 	q 			> 	0
			return p,q
		end
	end
	return 0,0
end

function firstbcdiff(A,B;maxdim=1,offset=0) # stands for first barcode difference
	for	r 			= 		0:maxdim
		Ba 			= 		sortslices(barcode(A,dim=r),		dims=1)
		Bb 			= 		sortslices(barcode(B,dim=r+offset),	dims=1)
		if 	Ba 		!= 		Bb
			return 	r
		end
	end
	return 0
end

function generatorbdc(C;dim=0)
	# bdc stands for birth, death, cycle status
	passedtest 		= 	true
	maxdim 			= 	C["input"]["maxdim"]
	B 				= 	barcode(C,dim=dim)
	if C["input"]["record"] == 	"all"
		for p 		= 	1:size(B,1)
			rep 	= 	classrep(C,dim=dim,class=p,format="index")
			case1 	= 	birthtime(C,dim=dim,chain=rep) == B[p,1]
			case2 	= 	deathtime(C,dim=dim,chain=rep) == B[p,2]
			case3 	= 	isempty(chainboundary(C,dim=dim,chain=rep))
			if 		! 	(case1 & case2 & case3)
				println([case1 case2 case3])
				println(rep)
				println([birthtime(C,dim=dim,chain=rep) B[p,1]])
				println([deathtime(C,dim=dim,chain=rep) B[p,2]])
				passedtest 	= 	false
				break
			end
		end
	else
		println()
		println("error in function <generatorbdc>: deathtimes can only be checked when C[\"input\"][\"record\"] = \"all\".")
		return
	end
	return passedtest
end

function ceilvr(	N;
					minrad 			= 	-Inf,
					maxrad 			= 	Inf,
					numrad 			= 	Inf)

	if minrad 						== 	"minedge"
		minrad 						= 	minimum(offdiagmin(N))
	end

	minrad,maxrad,mingrid,maxgrid 	= 	ceil2grid_overflowparameters(
										N,
										minrad 	= 	minrad,
										maxrad 	= 	maxrad)
					S 				= 	ceil2grid_overflow(
										N,
										minrad 	= 	minrad,
										maxrad 	= 	maxrad,
										mingrid = 	mingrid,
										maxgrid = 	maxgrid,
										numgrid = 	numrad)
	return S
end
