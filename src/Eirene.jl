# WELCOME TO EIRENE!
#
# You should have received a copy of the GNU General Public License along with
# Eirene.  If not, please see <http://www.gnu.org/licenses/>.
#
# Eirene Library for Homological Algebra
# Copyright (C) 2016, 2017, 2018, 2019, 2020  Gregory Henselman
# www.gregoryhenselman.org
#
# Eirene is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Eirene is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Eirene.  If not, see <http://www.gnu.org/licenses/>.
#
# PLEASE HELP US DOCUMENT Eirene's recent work! Bibtex entries and
# contact information for teaching and outreach can be found at the
# Eirene homepage, http://gregoryhenselman.org/eirene.

__precompile__()

module Eirene

##########################################################################################

#### 	REQUIREMENTS

##########################################################################################

using Pkg
using Distances
using JLD
using Blink
using PlotlyJS
using MultivariateStats
using Colors
using SparseArrays
using LinearAlgebra
using Dates
using Statistics
using DelimitedFiles
using CSV
using Hungarian #added for the Wasserstein distances

##########################################################################################

#### 	USER TEST FUNCTION

##########################################################################################

function example_function()
	print("Welcome to Eirene!  Great job running the example function.")
end


##########################################################################################

#### 	EXPORTS

##########################################################################################

export 	eirene,
		eirenefilepath,
		ezread,
		ezplot_pjs,
		plotpersistencediagram_pjs,
		plotclassrep_pjs,
		plotbarcode_pjs,
		plotbetticurve_pjs,
		ezplot_pjs,
		barcode,
		betticurve,
		classrep,
		latlon2euc,
		eirenefilepath,
		noisycircle,
		noisycircle3,
		torus,
		noisytorus,
		sphere,
		matchingcomplex_symmat,
		chessboardcomplex_symmat,
		plane2torus,
		zerodrandmat,
		ezlabel,
		unittest,
		wasserstein_distance #this is in "wassterstein_distances.jl"

##########################################################################################

#### 	SIMPLICIAL CONSTRUCTIONS

##########################################################################################

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0
	loci::Array{Int64,1} = copy(facenames)
	vrealization = Array{Int64}(undef,facecardinality,numfaces)
	post0::Int64 = 1
	post1::Int64 = 1

	for sd = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[sd]
		for i = 1:numfaces
			locus = loci[i]
			if cp[post0] > locus
				while cp[post0] > locus
					post0-=1
				end
				post1 = post0+1
			elseif cp[post1] <= locus
				while cp[post1] <= locus
					post1+=1
				end
				post0 = post1-1
			end
			loci[i] = farfaces[sd][locus]
			vrealization[sd,i]=post0
		end
	end
	return vrealization
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict,facecardinality,facenames)
	return vertexrealization(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict;dim = 1, class = 1)
	sd = dim+2
	facecard = dim+1

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		cyclename = barname2cyclename(D,class;dim = dim)
		rep = getcycle(D,sd,cyclename)
	end

	vrealization = vertexrealization(D::Dict,facecard,rep)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(farfaces,firstv,facecardinality,facenames)
	numfaces::Int64 = length(facenames)
	fc::Int64 = facecardinality
	m::Int64 = length(firstv[2])-1
	preallocationspace = 0
	vsupp = falses(m)
	loci::Array{Int64,1} = copy(facenames)
	post0::Int64 = 1
	post1::Int64 = 1
	for sd = facecardinality:-1:1
		cp::Array{Int64,1} = firstv[sd]
		for i = 1:numfaces
			locus = loci[i]
			if cp[post0] > locus
				while cp[post0] > locus
					post0-=1
				end
				post1 = post0+1
			elseif cp[post1] <= locus
				while cp[post1] <= locus
					post1+=1
				end
				post0 = post1-1
			end
			loci[i] = farfaces[sd][locus]
			vsupp[post0]=true
		end
	end
	return findall(vsupp)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict,facecardinality,facenames)
	return incidentverts(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict;dim=1,class=1)
	facecardinality = dim+1

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][dim+2][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,facecardinality,class)
	end
	return incidentverts(D,facecardinality,rep)
end

function buildclosefromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr;facecard = size(lclosefaces,1)+1)
	m = length(hcolptr)-1
	n = length(hrowval)
	hclosefaces = Array{Int64}(undef,facecard,n)
	if n == 0
		return hclosefaces
	else
		rowdepth = facecard-1
		rosettacol = Array{Int64}(undef,maximum(lrowval))
		for i = 1:m
			rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
			for j = cran(hcolptr,i)
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,j]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[facecard,j] = rosettacol[lrowval[farface]]
			end
		end
		return hclosefaces
	end
end

function buildallfromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr,selectedcolumnindices;verbose=false)
	if verbose
		println("PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY")
	end
	m = length(hcolptr)-1
	numhigs = length(hrowval)
	numselected = length(selectedcolumnindices)
	rowdepth = size(lclosefaces,1)
	sd = rowdepth+1
	hclosefaces = Array{Int64}(undef,sd+1,numselected)
	if numselected == 0
		return hclosefaces
	end
	rosettacol = Array{Int64}(undef,maximum(lrowval))
	columnsupp = falses(numhigs)
	columnsupp[selectedcolumnindices].=true
	columnmarker = 0
	for i = 1:m
		rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
		for j = cran(hcolptr,i)
			if columnsupp[j]
				columnmarker+=1
				farface = hrowval[j]
				for k = 1:rowdepth
					hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
				end
				hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
				hclosefaces[sd+1,columnmarker] = farface
			end
		end
	end
	return hclosefaces
end

function buildclosefromfar(farfaces,firstv,sd)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(undef,sd,n)
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)]=i
	end
	if sd == 2
		return lclosefaces'
	end
	for i = 3:sd
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	return lclosefaces
end

function buildclosefromfar(farfaces,firstv,sd,columnsinorder)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)].=i
	end
	if sd == 2
		return lclosefaces[columnsinorder]'
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
	end
	lclosefaces = buildclosefromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;facecard = sd-1)
	return lclosefaces
end

function buildallfromfar(farfaces,firstv,sd,columnsinorder;verbose = false)
	m = length(firstv[1])-1
	n = length(farfaces[sd])
	# destinationmatrix = Array{Int64}(undef,sd,n)
	if sd == 1
		return Array{Int64}(undef,0,m)
	end
	lclosefaces = Array{Int64}(undef,1,firstv[2][end]-1)
	for i = 1:m
		lclosefaces[cran(firstv[2],i)].=i
	end
	if sd == 2
		return vcat(lclosefaces[columnsinorder]',farfaces[sd][columnsinorder]')
	end
	for i = 3:(sd-1)
		lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
		#gc()
	end
	lclosefaces = buildallfromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;verbose = verbose)
	#gc()
	return lclosefaces
end

function ff2boundary(farfaces,firstv;sd=1)
	rv = Array{Int64}(undef,0)
	cp = [1]
	if sd == 1
		rv = Array{Int64}(undef,0)
		cp = ones(Int64,length(farfaces[1])+1)
	else
		n = length(farfaces[sd])
		rv = ff2aflight(farfaces,firstv,sd,1:n)
		rv = vec(rv)
		cp = convert(Array{Int64,1},sd+1:sd:(1+n*sd))
		prepend!(cp,[1])
	end
	return rv,cp
end

function ff2complex(farfaces,firstv;maxsd = length(farfaces))
	Nrv 	= fill(Array{Int64}(undef,0),maxsd)
	Ncp 	= fill(Array{Int64}(undef,0),maxsd)
	Nrv		= convert(Array{Array{Int64,1}},Nrv)
	Ncp		= convert(Array{Array{Int64,1}},Ncp)
	Nrv[1] 	= Array{Int64}(undef,0)
	Ncp[1]	= fill(1,length(farfaces[1])+1)
	for sd = 2:maxsd
		Nrv[sd],Ncp[sd] = ff2boundary(farfaces,firstv,sd=sd)
	end
	return Nrv,Ncp
end

function eirened2complex(C)
	if in(C["input"]["model"],["pc","vr"])
		rv,cp 	= 	boundarymatrices(C)
		fv 		= 	ocff2of(C["grain"],C["ocg2rad"])
	elseif in(C["input"]["model"],["complex"])
		rv 		= 	C["rv"]
		cp 		= 	C["cp"]
	else
		println("Error: the value of C[\"input\"][\"model\"] must be \"pc\", \"vr\", or \"complex\".")
	end
	fv 			= 	ocff2of(C["grain"],C["ocg2rad"])
	return 		rv,cp,fv
end

function ocff2of(grain::Array{Int64},ocg2rad::Array{Int64})
	m = length(grain)
	filt = Array{Int64}(undef,m)
	for i = 1:m
		filt[i] = ocg2rad[grain[i]]
	end
	return filt
end

function ff2aflight_sc2(farfaces,firstv,columns)
	sd = 2
	if isempty(farfaces[sd])
		return Array{Int64}(undef,2,0)
	end
	f0faces::Array{Int64,1} = farfaces[sd]
	colptr::Array{Int64,1} = firstv[2]
	columnpost::Int64   = 1
	columnpostp1::Int64 = 2
	faces::Array{Int64,2} = Array{Int64}(undef,2,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		if f0 >= colptr[columnpostp1]
			while f0 >= colptr[columnpostp1]
				columnpostp1+=1
			end
			columnpost = columnpostp1-1
		elseif f0 < colptr[columnpost]
			while f0 < colptr[columnpost]
				columnpost-=1
			end
			columnpostp1 = columnpost+1
		end
		faces[1,fp] = columnpost
		faces[2,fp] = f0faces[f0]
	end
	return faces
end

function ff2aflight_sc3(farfaces,firstv,columns)
	sd = 3

	if isempty(farfaces[sd])
		return Array{Int64}(undef,3,0)
	end

	fcfaces::Array{Int64,2} = buildclosefromfar(farfaces,firstv,sd-1,1:length(farfaces[2]))

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]

	holdi=[1];holdip1=[2]
	t1::Array{Int64,1} = Array{Int64}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)

	faces::Array{Int64,2} = Array{Int64}(undef,3,length(columns))
	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		f3 = fcfaces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1} ,holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		faces[1,fp] = t1[f3]
		faces[2,fp] = t1[f2]
		faces[3,fp] = f1
	end
	return faces
end

function ff2aflight_scgt3(farfaces,firstv,sd,columns)

	if isempty(farfaces[sd])
		return Array{Int64}(undef,sd,0)
	end

	f0faces::Array{Int64,1} = farfaces[sd]
	f1faces::Array{Int64,1} = farfaces[sd-1]
	f2faces::Array{Int64,1} = farfaces[sd-2]
	fcfaces::Array{Int64,2} = buildallfromfar(farfaces,firstv,sd-2,1:(firstv[sd-2][end]-1))

	fvscm0::Array{Int64,1}  = firstv[sd]
	fvscm1::Array{Int64,1}  = firstv[sd-1]
	fvscm2::Array{Int64,1}  = firstv[sd-2]
	fvscm3::Array{Int64,1}  = firstv[sd-3]

	holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
	t1::Array{Int64,1} = Array{Int64}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
	t2::Array{Int64,1} = Array{Int64}(undef,fvscm3[end]-1);t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

	scm0::Int64 = sd; scm1::Int64 = sd-1; scm2::Int64 = sd-2
	faces::Array{Int64,2} = Array{Int64}(undef,sd,length(columns))

	for fp = 1:length(columns)
		f0 = columns[fp]
		f1 = f0faces[f0]
		f2 = f1faces[f1]
		updatetranslator!(f0::Int64,fvscm0::Array{Int64,1},holdi::Array{Int64,1},holdip1::Array{Int64,1},t1::Array{Int64,1},fvscm1::Array{Int64,1},f1faces::Array{Int64,1})
		updatetranslator!(f1::Int64,fvscm1::Array{Int64,1},holdj::Array{Int64,1},holdjp1::Array{Int64,1},t2::Array{Int64,1},fvscm2::Array{Int64,1},f2faces::Array{Int64,1})
		for i = 1:scm2
			faces[i,fp] = t1[t2[fcfaces[i,f2]]]
		end
		faces[scm1,fp] = t1[f2]
		faces[scm0,fp] = f1
	end
	return faces
end

# hdlip1 is always used with value 2- maybe better make it local or with default value?
function updatetranslator!(f0,firstv0,holdi,holdip1,t,firstv1,farfaces1)
	if firstv0[holdip1[1]] <= f0
		while firstv0[holdip1[1]]<= f0
			holdip1[1]+=1
		end
		holdi[1] = holdip1[1]-1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	elseif firstv0[holdi[1]] > f0
		while firstv0[holdi[1]] > f0
			holdi[1]-=1
		end
		holdip1[1] = holdi[1]+1
		t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
	end
end


function ff2aflight(farfaces,firstv,sd,columns)
	if sd == 1
		return Array{Int64}(undef,0,length(columns))
	elseif sd == 2
		return ff2aflight_sc2(farfaces,firstv,columns)
	elseif sd == 3
		return ff2aflight_sc3(farfaces,firstv,columns)
	else
		return ff2aflight_scgt3(farfaces,firstv,sd,columns)
	end
end

function ff2aflight(D::Dict,sd,columns)
	farfaces = D["farfaces"]; firstv = D["firstv"]
	faces = ff2aflight(farfaces,firstv,sd,columns)
	return faces
end

#=

NB
- Input argument <grain> must be arranged least to greatest

OUTPUTS
- higlab
The concatenated vector [pphigs,nphigs)]
- lowlab
The concatenated vector [pplows,nplows[perm]], where perm is a permutation such
that the entries of lowgrain[nplows[perm]] appear in ascending order, numer-
ically.
- Mrv, Mcp, Mm
Sparse matrix representation of transpose(D[lowlab,higlab]), where D is
submatrix of the total boundary operator indexed by cells of dimension sd-1
(along the columns) and sd-2 (along the rows).

=#
function filteredmatrixfromfarfaces(
	farfaces,
	firstv,
	prepairs,
	grain,
	sd::Integer,
	lowbasisnames::Array{Int64,1};
	verbose = false)

	numhigs = length(farfaces[sd])
	numlows = length(farfaces[sd-1])
	numppair= length(prepairs[sd])

	pphigs = prepairs[sd]
	pplows = farfaces[sd][pphigs]
  	lpls = lowbasisnames
	hphs = farfaces[sd+1][prepairs[sd+1]]
	nplows = intervalcomplementuniqueunsortedinput(vcat(lpls,pplows),numlows)
	nphigs = intervalcomplementuniqueunsortedinput(vcat(hphs,pphigs),numhigs)

	numnhph = numhigs-length(hphs)
	Ml = numlows - length(lpls)
	Mh = numhigs - length(hphs)

	higtranslator = zeros(Int64,numnhph)
	lowtranslator = zeros(Int64,numlows)
	lowtranslator[pplows] = 1:numppair

	# if !isempty(nplows) && sd > 2
	# 	npfilt = grain[sd-1][nplows]
	# 	nporder = integersinsameorder(npfilt)
	# 	addinteger!(nporder,numppair)
	# else
	# 	nporder = (numppair+1):(numppair+length(nplows))
	# end

	if sd > 1
		npfilt 	= 	grain[sd-1][nplows]
		nporder = 	integersinsameorder(npfilt)
		addinteger!(nporder,numppair)
	else
		npfilt 	= 	zeros(Int64,0)
		nporder =	zeros(Int64,0)
	end

	lowtranslator[nplows] = nporder
	higsinpointorder = intervalcomplementuniqueunsortedinput(hphs,numhigs)
	lowlab = Array{Int64}(undef,Ml)
	lowlab[1:numppair]=pplows
	lowlab[nporder]=nplows
	higlab = vcat(pphigs,nphigs)

	if verbose
		comparisonsuppvec = trues(numhigs)
		comparisonsuppvec[hphs]=false
		comparisonvec=findall(comparisonsuppvec)
		differencecounter = 0
		for i = 1:length(higsinpointorder)
			if higsinpointorder[i]!=comparisonvec[i]
				differencecounter+=1
			end
		end
		if differencecounter>0
			print(["hi ho comparison vec" differencecounter])
			print(length(higsinpointorder))
			print(length(comparisonvec))
			print(comparisonvec[1:20])
			print(higsinpointorder[1:20])
			sleep(5)
		end
	end
	ppsupp = falses(numhigs)
	ppsupp[pphigs].=true
	ppmarker = 0
	nppmarker = numppair
	for i = 1:numnhph
		hig = higsinpointorder[i]
		if ppsupp[hig]
			ppmarker+=1
			higtranslator[i]=ppmarker
		else
			nppmarker+=1
			higtranslator[i]=nppmarker
		end
	end
	allfaces = buildallfromfar(farfaces,firstv,sd,higsinpointorder;verbose = verbose)
	if verbose
		print("done building allfromfar")
	end
	Mrv,Mcp,Mm = presparsefull2unsortedsparsetranspose(allfaces,lowtranslator,higtranslator;verbose=verbose)
	higtranslator = [];npfilt = [];ppsupp = [];allfaces = []
	#gc()
	if verbose && length(Mrv)>(Mcp[end]-1)
		print("There was the thought that Mrv should have no extra elements")
		sleep(3)
	end
	return Mrv,Mcp,lowlab,higlab,Mm
end


function grain2maxsd(grain)
	c = 0
	for i = 1:length(grain)
		if !isempty(grain[i])
			c = i
		end
	end
	return c
end

##########################################################################################

####	CELL OPERATIONS

##########################################################################################

function ocff2of(grain::Array{Int64},ocg2rad::Array{Float64})
	m = length(grain)
	filt = Array{Float64}(undef,m)
	for i = 1:m
		filt[i] = ocg2rad[grain[i]]
	end
	return filt
end

function ocff2of(grain::Array{Array{Int64,1},1},ocg2rad::Array{Float64})
	n = length(grain)
	filt = Array{Array{Float64}}(undef,n)
	for i = 1:n
		filt[i] = ocff2of(grain[i],ocg2rad)
	end
	return filt
end
##########################################################################################

####	SCHUR COMPLEMENTS

##########################################################################################

function schurit4!(	Mrv,Mcp,Mm,Mn,Mn0,
						rowlab,collab,
						Jprows,Jpcols,numjunpairs,
						Sprows,Spcols,numsenpairs,
						comprows,compcols,
						Trv,Tcp,Srv,Scp;
						updatetransform = true,
						verbose = false,
						diagonstic = false
					)

	Mm0 = copy(Mm[1])
	Mm[1] = length(comprows)
	Mn[1] = length(compcols)

	copycolind2colind!(Trv,Tcp,Jpcols[numjunpairs[1]:-1:1],Srv,Scp,numsenpairs[1]+1,0)

	topspot = numsenpairs[1]+numjunpairs[1]
	for i = 1:numjunpairs[1]
		placementlocation = numsenpairs[1]+i
		extractionlocation = numjunpairs[1]-i+1
		Sprows[placementlocation] = rowlab[Jprows[extractionlocation]]
		Spcols[placementlocation] = collab[Jpcols[extractionlocation]]
	end

	numsenpairs[1]+=numjunpairs[1]

	rowsum = falses(Mm0)#zeros(Int64,Mm0)
	for jp = 1:Mn[1]
		j = compcols[jp]
		for ip = cran(Mcp,j)
			rowsum[Mrv[ip]]=true
		end
	end

	keptlist = finddownstreamelements_embeddedupperunitriangularmatrix(
					Mrv,Mcp,Mm0,findall(rowsum),Jprows[1:numjunpairs[1]],Jpcols[1:numjunpairs[1]];verbose=verbose
					)

	if verbose
		println()
		println([length(keptlist) "=numkept" numjunpairs[1] "=numinputp" Mn[1] "=Mn" Mm[1] "=Mm" (Mcp[Mn[1]+1]-1) "=nnz(M)"])
	end

	keptmarker = length(keptlist)
	prows = Array{Int64}(undef,keptmarker)
	pcols = Array{Int64}(undef,keptmarker)
	for i = 1:keptmarker
		keptindex = keptlist[i]
		prows[i] = Jprows[keptindex]
		pcols[i] = Jpcols[keptindex]
	end

	Arv,Crv,Acp,Ccp = stackedsubmatrices(Mrv,Mcp,prows,comprows,pcols,Mm0)
	Brv,Drv,Bcp,Dcp = stackedsubmatrices(Mrv,Mcp,prows,comprows,compcols,Mm0)
	Lrv,Lcp = copycolumnsubmatrix(Trv,Tcp,pcols)
	Rrv,Rcp = copycolumnsubmatrix(Trv,Tcp,compcols)

	translator = Array{Int64}(undef,Mm0)
	translator[prows]=1:keptmarker
	yafterx!(translator,Arv)
	yafterx!(translator,Brv)

	translator[comprows]=1:length(comprows)
	yafterx!(translator,Crv)
	yafterx!(translator,Drv)

	Airv,Aicp = morseInverseF2orderedColsUnsortedRowsInSilentOut(Arv,Acp)
	Brv,Bcp = spmmF2silentLeft(Airv,Aicp,Brv,Bcp,keptmarker)

	for j = 1:keptmarker
		translator[j]=collab[pcols[j]] #repurposing rowsum as col-to-row translator
	end
	collabcopy = copy(collab)
	for i = 1:Mm[1]
		rowlab[i] = rowlab[comprows[i]]
	end
	for j = 1:Mn[1]
		collab[j] = collab[compcols[j]]
	end

	blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Mm[1],Mn[1],Mrv,Mcp,0)

	if updatetransform
		blockprodsumsilenticolsleftWrite2!(Lrv,Lcp,Brv,Bcp,Rrv,Rcp,Mn0,Mn[1],Trv,Tcp,0,translator)
	end
end

##########################################################################################

####	CHAIN OPERATIONS

##########################################################################################

#=
notes on morselu!
- 	the output array tid has size equal to the number of columns of the
	input array, NOT the column-rank of the input array
-	the columns of M must be ordered by grain (ascending)
-	the first (rank of M) elements of tlab index the complete set
	of nonzero columns in the reduced matrix
=#
function morselu!(
	Mrv::Array{Tv,1},
	Mrowgrain::Array{Tv,1},
	Mcp::Array{Tv,1},
	Mcolgrain::Array{Tv,1},
	lowlab::Array{Tv,1},
	higlab::Array{Tv,1},
	pplow::Array{Tv,1},
	pphig::Array{Tv,1},
	Mm::Integer;
	storetransform = true,
	verbose = false,
	diagnostic = false) where Tv<:Integer

 	rowlab = higlab;collab = lowlab

	Mm = [length(higlab)]
	Mn = [length(lowlab)]
	Mn0 = Mn[1]
	maxnz = Mcp[Mn[1]+1]

	maxnumpairs = min(Mm[1],Mn[1]); numjunpairs = [length(pplow)]; numsenpairs = [0]
	Sprows=Array{Tv}(undef,maxnumpairs);Spcols=Array{Tv}(undef,maxnumpairs);
	Jprows=Array{Tv}(undef,maxnumpairs);Jpcols=Array{Tv}(undef,maxnumpairs);
	Jprows[1:numjunpairs[1]]=pphig;Jpcols[1:numjunpairs[1]]=pplow
	comprows = convert(Array{Tv,1},(numjunpairs[1]+1):Mm[1])
	compcols = convert(Array{Tv,1},(numjunpairs[1]+1):Mn[1])

	Trv=Array{Tv}(undef,0);Srv = Array{Tv}(undef,0)
	Tcp=ones(Tv,Mn[1]+1);Scp=ones(Tv,Mn[1]+1)

	if diagnostic
		numsenpairsOLD = numsenpairs[1]
	end
	schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
					rowlab,collab,
					Jprows,Jpcols,numjunpairs,
					Sprows,Spcols,numsenpairs,
					comprows,compcols,
					Trv,Tcp,Srv,Scp;
					updatetransform = storetransform,
					verbose = verbose
				)
	maxnz = max(maxnz,Mcp[Mn[1]+1])

	if verbose
		println("first shurr finished")
		if Mn[1]>0
			println([Mcp[Mn[1]] "nnz(M)" Mm[1] "Mm" Mn[1] "Mn" length(Mrv) "length(Mrv)"])
		else
			println("Mn = 0")
		end
	end
	#gc()
	rowfilt = Array{Tv}(undef,length(comprows)); colfilt = Array{Tv}(undef,length(compcols))
	counter = 0
	while Mcp[Mn[1]+1]>1
		if verbose
			println("starting block operation $(counter)")
		end
		counter+=1
		for i = 1:Mm[1]
			rowfilt[i] = Mrowgrain[rowlab[i]]
		end
		for j = 1:Mn[1]
			colfilt[j] = Mcolgrain[collab[j]]
		end
		getPairsLightWrite2!(Mrv,Mcp,rowfilt,colfilt,Mm[1],Mn[1],Jprows,Jpcols,numjunpairs,verbose=verbose)
		comprows = intervalcomplementuniqueunsortedinput(Jprows[1:numjunpairs[1]],Mm[1])
		compcols = intervalcomplementuniqueunsortedinput(Jpcols[1:numjunpairs[1]],Mn[1])
		if diagnostic
			numsenpairsOLD = numsenpairs[1]
		end

		schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
						rowlab,collab,
						Jprows,Jpcols,numjunpairs,
						Sprows,Spcols,numsenpairs,
						comprows,compcols,
						Trv,Tcp,Srv,Scp;
						updatetransform = storetransform,
						verbose = verbose
					)
		maxnz = max(maxnz,Mcp[Mn[1]+1])
	end
	lastSrowmarker = Scp[numsenpairs[1]+1]-1
	lastTrowmarker = Tcp[Mn[1]+1]-1
	deleteat!(Srv,(lastSrowmarker+1):length(Srv))
	deleteat!(Trv,(lastTrowmarker+1):length(Trv))
	deleteat!(Scp,(numsenpairs[1]+1):length(Scp))
	deleteat!(Tcp,(Mn[1]+2):length(Tcp))
	deleteat!(Sprows,(numsenpairs[1]+1):maxnumpairs)
	deleteat!(Spcols,(numsenpairs[1]+1):maxnumpairs)
	Tcp.+=lastSrowmarker
	append!(Scp,Tcp)
	append!(Srv,Trv[1:lastTrowmarker])
  	tlab = Spcols[1:numsenpairs[1]]
	append!(tlab,collab[1:Mn[1]])
	return Srv,Scp,Sprows,Spcols,tlab,maxnz
end

function persistf2_core_cell(
	Nrv,
	Ncp,
	grain;
	maxsd = length(Nrv),
	record="cyclerep",
	verbose=false,
	prepairs = fill(Array{Int64}(undef,0),maxsd+1)
	)
	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1)
	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1)
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1)
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1)
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1)
	for i in [1,maxsd+1]
		tcp[i] 				=[1]
		trv[i]				=Array{Int64}(undef,0)
		phi[i]				=Array{Int64}(undef,0)
		plo[i]				=Array{Int64}(undef,0)
		tid[i]				=Array{Int64}(undef,0)
	end
	maxnzs 									=zeros(Int64,maxsd+1)
	m 										=length(Ncp[1])-1
	for sd = 2:maxsd
		if sd > length(Nrv)
			trv[sd] = Array{Int64}(undef,0)
			tcp[sd] = ones(Int64,1)
			tid[sd] = Array{Int64}(undef,0)
			plo[sd] = Array{Int64}(undef,0)
			phi[sd] = Array{Int64}(undef,0)
			continue
		elseif sd>2
			lowbasisnames = phi[sd-1]
		else
			lowbasisnames = Array{Int64}(undef,0)
		end
		Mrv 			= Nrv[sd] # temporary
		Mcp 			= Ncp[sd] # temporary
		if isempty(Mrv)
			plo[sd] 	= Array{Int64}(undef,0)
			phi[sd] 	= Array{Int64}(undef,0)
			tid[sd] 	= Array{Int64}(1:numcols(Ncp[sd-1]))
						  deleteat!(tid[sd],sort(phi[sd-1]))
		    perm 		= sortperm(grain[sd-1][tid[sd]],alg=MergeSort)
			tid[sd] 	= tid[sd][perm]
			trv[sd] 	= Array{Int64}(undef,0)
			tcp[sd] 	= ones(Int64,1+length(tid[sd]))
			continue
		end
		Mm0				= length(Ncp[sd-1])-1
		Mn0 			= length(Ncp[sd  ])-1
		Mm				= [Mm0]  					# temporary
		Mn 				= [Mn0] 					# temporary
		higlab			= convert(Array{Int64},1:Mn[1])	# we'll assume prepairs to be empty, for now
		lowlab	 		= intervalcomplementuniqueunsortedinput(lowbasisnames,Mm0)	# temporary, and we'll assume prepairs is empty for now
		nporder			= sortperm(grain[sd-1][lowlab],alg=MergeSort)
		lowlab			= lowlab[nporder]
		Mrv,Mcp			= transposeLighter_submatrix(
							Mrv, 					# the Arv argument
							Mcp, 					# the Acp argument
							Mm0, 					# the Am argument
							rows = lowlab,			# the rows selected
							cols = higlab)			# the columns selected
 		lowlabtemp 		= convert(Array{Int64,1},1:length(lowlab))
 		higlabtemp 		= convert(Array{Int64,1},1:length(higlab))
 		higfilttemp 	= grain[sd][higlab]
 		lowfilttemp 	= grain[sd-1][lowlab]
		pplow 			= convert(Array,length(prepairs[sd]):-1:1) # we'll assume prepairs is empty, for now
		pphig 			= convert(Array,length(prepairs[sd]):-1:1) # we'll assume prepairs is empty, for now
		if verbose
			println("Constructed Morse boundary operator, columns indexed by cells of dimension $(sd-1)")
		end
		# NB: It is critical that the columns of the input array should
		# be ordered according to filtration; in particular, the entries of
		# lowfiltemp should increase monotonically
		Srv,Scp,Sphigs,Splows,tlab,maxnz =
		morselu!(
			Mrv,
			higfilttemp,
			Mcp,
			lowfilttemp,
			lowlabtemp,
			higlabtemp,
			pplow,
			pphig,
			Mm[1],
			storetransform = storetransform,
			verbose = verbose)
		trv[sd] 		= lowlab[Srv]
		tcp[sd] 		= Scp
		tid[sd] 		= lowlab[tlab]
		plo[sd] 		= lowlab[Splows]
		phi[sd] 		= higlab[Sphigs]
		maxnzs[sd] 	= maxnz
	end
	return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2_core_vr(
	farfaces::Array{Array{Int64,1},1},
	firstv::Array{Array{Int64,1},1},
	prepairs::Array{Array{Int64,1},1},
	grain::Array{Array{Int64,1},1},
	maxsd::Integer;
	record="cyclerep",
	verbose=false)

	farfaces::Array{Array{Int64,1},1}
	firstv::Array{Array{Int64,1},1}
	prepairs::Array{Array{Int64,1},1}
	grain::Array{Array{Int64,1},1}

	if record == "all" || record == "cyclerep"
		storetransform = true
	else
		storetransform = false
	end

	m = length(firstv[1])-1

	trv::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1);
	tcp::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1);
	phi::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1);
	plo::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1);
	tid::Array{Array{Int64,1},1}			=Array{Array{Int64,1},1}(undef,maxsd+1);
	for i in [1,maxsd+1]
		tcp[i]  			=[1];
		trv[i]				=Array{Int64}(undef,0)
		tid[i]				=Array{Int64}(undef,0)
		phi[i]				=Array{Int64}(undef,0)
		plo[i]				=Array{Int64}(undef,0)
	end

	maxnzs 									=zeros(Int64,maxsd+1);

	for sd = 2:maxsd
		if sd > length(farfaces)
			continue
		elseif sd>2
			lowbasisnames = phi[sd-1]
		else
			lowbasisnames = Array{Int64}(undef,0)
		end
		Mrv::Array{Int64,1},
		Mcp::Array{Int64,1},
		lowlab::Array{Int64,1},
		higlab::Array{Int64,1},
		Mm =
		filteredmatrixfromfarfaces(farfaces,firstv,prepairs,grain,sd,lowbasisnames;verbose=verbose)
 		lowlabtemp = convert(Array{Int64,1},1:length(lowlab))
 		higlabtemp = convert(Array{Int64,1},1:length(higlab))
 		higfilttemp = grain[sd][higlab]
 		lowfilttemp = grain[sd-1][lowlab]
		if verbose
			println("Constructed Morse boundary operator, columns indexed by cells of dimension $(sd-1)")
		end
		pplow = convert(Array,length(prepairs[sd]):-1:1)
		pphig = convert(Array,length(prepairs[sd]):-1:1)

		# NB: It is critical that the columns of the input array
		# be ordered according to filtration; in particular, the entries of
		# lowfiltemp should increase monotonically
		Srv,Scp,Sphigs,Splows,tlab,maxnz =
		morselu!(
			Mrv,
			higfilttemp,
			Mcp,
			lowfilttemp,
			lowlabtemp,
			higlabtemp,
			pplow,
			pphig,
			Mm,
			storetransform = storetransform,
			verbose = verbose)
		trv[sd] = lowlab[Srv]
		tcp[sd] = Scp
		tid[sd] = lowlab[tlab]
		plo[sd] = lowlab[Splows]
		phi[sd] = higlab[Sphigs]
		maxnzs[sd]= maxnz
	end
	return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2!(
	D::Dict;maxsd=0,
	dictionaryoutput::Bool = true,
	verbose::Bool = false,
	record = "cyclerep")

	farfaces = D["farfaces"]
	firstv = D["firstv"]
	prepairs = D["prepairs"]
	grain = D["grain"]
	if maxsd == 0
		maxsd = length(farfaces)-1
	end

	trv,tcp,plo,phi,tid,maxnzs =
	persistf2_core_vr(farfaces,firstv,prepairs,grain,maxsd::Integer;record=record,verbose = verbose)
	if dictionaryoutput == true
		D["trv"] = trv
		D["tcp"] = tcp
		D["tid"] = tid
		D["plo"] = plo
		D["phi"] = phi
		D["maxnz"] = maxnzs
		return D
	else
		return trv,tcp,plo,phi,tid
	end
end

function persistf2vr(
	s,
	maxsd;
	model 			= "vr",
	entryformat 		= "textfile",
	minrad			= -Inf,
	maxrad			= Inf,
	numrad			= Inf,
	nodrad  		= [],
	filfun 			= "n/a", # stands for filtration function; only to be used with fastop=false; function must take finite values; not compatible with minrad, maxrad, or numrad
	fastop			= true,
	vscale			= "diagonal",
	pointlabels 	= [],
	verbose 		= false,
	record 			= "cyclerep")

	#### Start timer
# 	tic()

	#### Must not stop early if an explicit filtration function is passed; also do no rounding
	if filfun 	   != 	"n/a"
		fastop 		= 	false
		numrad 		= 	1
		println("must define <clearprepairs>")
		return
	end

	#### Extract data as necessary
	inputisfile = false
	if typeof(s) == String
		inputisfile = true
		filename = s #modified 12/29/2017; note that strings are immutable
		if entryformat == "textfile"
			if typeof(readdlm(filename,','))<:Array{Float64}
				s = readdlm(s,',')
			elseif typeof(readdlm(filename,' '))<:Array{Float64}
				s = readdlm(s,' ')
			else
				print("Error reading text file.  Input files of this format must be either comma or space delimited.")
				return
			end
		elseif entryformat == "perseus_distmat"
			s,minrad,maxrad,maxsd = parse_perseusdistmat(filename)
			model = "clique_perseus"
			numrad = Inf
		elseif entryformat == "perseus_brips"
			s,minrad,maxrad = parse_perseusbrips(filename)
			s = s'  # the perseus brips format stores points as rows
			model = "pointcloud_perseusbrips"
			numrad = Inf
		end
	else
		filename = "user input julia array"
	end
	if model == "pc" || model == "perseus_brips"
		pc 	= 	"genera"
	else
		pc 	= 	"n/a"
	end

	#### Store the input
	input = Dict(
		"model"			=> model,
		"genera"		=> copy(s),
		"pc"			=> pc,
		"source"		=> filename,
		"maxdim" 		=> maxsd-2,
		"maxrad"		=> maxrad,
		"minrad"		=> minrad,
		"numrad"		=> numrad,
		"nodrad"		=> nodrad,
		"fastop"		=> fastop,
		"record"	 	=> record,
		"filfun" 		=> filfun,
		# "version" 		=> Pkg.installed("Eirene"),
		"date"			=> string(Dates.Date(now())),
		"time"			=> string(Dates.Time(now()))
		)

	#### Determine the number of points
	if model == "pc"
		numpoints = size(s,2)
	elseif model == "vr"
		numpoints = size(s,1)
		if !issymmetric(s)
			print("It appears the input matrix is not symmetric.  Only symmetric distance matrices are accepted when the <model> keyword argument has value \"clique\".")
			return
		end
	else
		print("Keyword argument <model> must be \"vr\", \"pc\", or \"cell\".  Please see documentaiton for further details.")
	end

	#### Extract labels
	if pointlabels == "none" || pointlabels == []
		pointlabels = 1:numpoints
	elseif pointlabels in ["left","right","top","bottom"]
		(s,pointlabels) = separatelabels(s,pointlabels)
	elseif typeof(pointlabels) == String
		input["source_pointlabels"] = pointlabels
		pointlabels = ezread(pointlabels)
	end
	if length(pointlabels) != numpoints
		warn("It appears the number of vertex labels does not match the number of vertices.")
	end
	pointlabels = ezlabel(pointlabels)
	for i = 1:length(pointlabels)
		pointlabels[i] = "$(pointlabels[i])"
	end
	input["pointlabels"] = pointlabels

	#### type the matrix
	s = convert(Array{Float64,2},s)
	if model == "pc"
		d = Distances.pairwise(Euclidean(),s,dims=2)
		if !isempty(nodrad)
			for 	i 	= 	1:numpoints
				d[i,i] 	= 	nodrad[i]
			end
		end
	elseif model == "vr"
		d = convert(Array{Float64,2},s)
	end

	################################################################################

	if fastop
		maxrad_alt 	= 	minimum(maximum(d,dims=1))
		maxrad_alt  = 	min(maxrad_alt,maxrad)
	else
		maxrad_alt 	= 	maxrad
	end

	d 			= 	minmaxceil(		d,
									minrad 	= 	minrad,
									maxrad 	= 	maxrad_alt,
									numrad 	= 	numrad)

	# vfilt 		= 	diag(t)
	# recall that t will have the same order (NOT inverted) as d
	# <trueordercanonicalform> is a bit like <integersinsameorder>, just valid for floating point inputs, and with a bit more data in the output
	t,ocg2rad 	= 	trueordercanonicalform(d,factor=true)

	t 			= 	(1+maximum(t)).-t
	ocg2rad 	= 	reverse(ocg2rad,dims=1)

	if 	any(d.>maxrad_alt)
		t 		= 	t.-1
		deleteat!(ocg2rad,1)
	end

	vertices2keep 	= 	findall(diag(t).!=0)  # this step is necessary in order to cover the case where some vertices never enter the filtration
	t 				= 	t[vertices2keep,vertices2keep]

	#### Build the complex
	D = buildcomplex3(t,maxsd;verbose = verbose)
	D["ocg2rad"]=ocg2rad

	################################################################################

	################################################################################
	# arbitrary function values
	# NB assumes the function takes only finite values
	if filfun != "n/a"

		fv 	= 	Array{Array{Float64,1}}(maxsd) # fv stands for filtration values
		for sd = 1:maxsd
			fv[sd] 	= 	Array{Float64}(undef,length(D["farfaces"][sd]))
			for p 	= 	1:length(D["grain"][sd])
				# syntax reminder: vertexrealization(D::Dict,facecardinality,facenames)
				fv[sd][p] 	= 	filfun(vertexrealization(D,sd,[p]))
			end
		end

		ocg,ocg2rad 	= 	trueordercanonicalform(cat(fv...,dims=1),factor=true) # ocg stands for order canonical grain
		ocg 			= 	1+maximum(ocg)-ocg
		ocg2rad 		= 	reverse(ocg2rad,dims=1)


		D["ocg2rad"]	=	ocg2rad
		D["grain"] 		= 	ocg
	end

	################################################################################

	#### Compute persistence
	persistf2!(D;verbose = verbose,record = record)

	#### Store input data
	D["input"] 		= 	input
	D["nvl2ovl"]	= 	vertices2keep[D["nvl2ovl"]]  # this covers the case where some vertices never enter the filtration

	#### Store generators
	#gc()
	if record == "all" || record == "cyclerep"
 		unpack!(D)
	end
	#gc()
	if record == "cyclerep"
		delete!(D,"trv")
		delete!(D,"tcp")
		delete!(D,"Lrv")
		delete!(D,"Lcp")
		delete!(D,"Rrv")
		delete!(D,"Rcp")
		delete!(D,"Lirv")
		delete!(D,"Licp")
		delete!(D,"prepairs")
	end

	#### Record time
# 	D["input"]["computationtime"] = toc()

	return D
end

#=
- 	drafted 12/27/2017
- 	version with 1 (string) non-keyword argument
- 	formatting guidelines:
	- the input file should be a csv in which
	- line 1: cell dimensions
	- line 2: cell filtrations
	- line 3: boundary matrix row values
	- line 4: boundary matrix column pattern

NB: the default value for maxdim has not been tested, and may cause errors;
the -3 accounts for (1) julia uses 1-indexed arrays, (2) to calculate
homology in dimension p, one must inspect the (p+1)-dimensional boundary
operator, (3) this operator should be givent the same treatment as those
that precede it ... generally this assumes that the next one up is at least
defined, even if it is trivial.
=#
function persistf2complex(filepath::String;
						maxdim=Inf,
						entryformat = "sp",
						record = "cyclerep",
						verbose=false)

	if entryformat 		== 	"sp"
		dp,fv,rv,cp 	= 	humanreadablefilepath2unsegmentedfilteredcomplex(filepath)
	elseif in(entryformat,["dp","dv","ev"])
		dp,fv,rv,cp 	= 	cscfilepath2unsegmentedfilteredcomplex(filepath,toprow = entryformat)
	end

	if 	maxdim 		== 	Inf
		maxdim 		= 	length(dp)-3
	end

	C  				= 	persistf2complex(
						rv = rv,
						cp = cp,
						fv = fv,
						dp = dp,
						maxdim=maxdim,
						record=record,
						verbose=verbose)
	C["input"]["source"]		= 	filepath;
	C["input"]["entryformat"]	= 	entryformat;
	return C
end

function emptyunsegmentedfilteredcomplex_dp()
	rv 	= zeros(Int64,0)
	cp 	= ones(Int64,1)
	fv	= zeros(Float64,0)
	dp  = ones(Int64,1)
	return dp,fv,rv,cp #rv,cp,fv,dp
end

# possible values for toprow: dp, dv, ev
function cscfilepath2unsegmentedfilteredcomplex(fp;toprow="dp")
	M 				= 	CSV.read(fp,header=0,silencewarnings=true)
	nemo 			=	Array{Any}(M[:,1])

	#	a zero operator // empty complex, formatted by [dv, ev]
	if isequal(nemo,[missing,missing,missing,1])
		if in(toprow,["dv","ev"])
			return emptyunsegmentedfilteredcomplex_dp() # this returns values for rv,cp,fv,dp
		else
			print("Error: please check formatting of input file.")
			return
		end

	#	a zero operator // empty complex, formatted by [dp]
	elseif isequal(nemo,[1,missing,missing,1])

		if in(toprow,["dp"])
			return emptyunsegmentedfilteredcomplex_dp() # this returns values for rv,cp,fv,dp
		else
			print("Error: input file should have .csv format with four lines.")
			return
		end

	#	a zero operator // positive # of cells, empty boundary
	elseif ismissing.(nemo) == [false,false,true,false]

		if any(M[3,:].!=1)
			print("Error: please check formatting of input file.")
			return
		else
			# dimension pattern
			xx 	    = csvreadrow(fp,row=1,rowtype=Int64)
			# filtration values
			fv      = csvreadrow(fp,row=1,rowtype=Float64)
			# row values
			rv 	    = zeros(Int64,0)
			# column pattern
			cp      = csvreadrow(fp,row=3,rowtype=Int64)
		end

	#	nonzero operator
	elseif ismissing.(nemo) == [false,false,false,false]
		# dimension pattern
		xx 	    	= 	csvreadrow(fp,row=1,rowtype=Int64)
		# filtration values
		fv      	= 	csvreadrow(fp,row=2,rowtype=Float64)
		# row values
		rv 	    	= 	csvreadrow(fp,row=3,rowtype=Int64)
		# column pattern
		cp      	= 	csvreadrow(fp,row=4,rowtype=Int64)

	#	error
	else
		println("Error: please check formatting of input file.")
		printval(xx,"missing values")

	end

	if toprow 			== 	"dp"
		dp 				= 	xx
	elseif toprow	 	== 	"dv"
		dp 				= 	dimensionvalues2dimensionpattern(xx)
	elseif toprow	 	== 	"ev"
		dp 				= 	eulervector2dimensionpattern(xx)
	end

	return dp,fv,rv,cp
end

function humanreadablefilepath2unsegmentedfilteredcomplex(fp)
	M 					= 	CSV.read(fp,header=0,silencewarnings=true)
	M 					= 	convert(Matrix{Float64},M)
	m 					= 	size(M,1)

	dv 					= 	Array{Int64}(M[:,1])
	fv 					= 	Array{Float64}(M[:,2])
	dp 					= 	Eirene.dimensionvalues2dimensionpattern(dv)

	rv 					=	zeros(Int64,0)
	cp 					= 	zeros(Int64,m+1)
	cp[1]				= 	1

	for p 				= 	1:m
		vals  			=	csvreadrow(fp,row=p,rowtype=Float64)
		vals 			=	vals[3:end]
		append!(rv,vals)
		cp[p+1]  		=	length(rv)+1
	end

	return dp,fv,rv,cp
end


function unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ncd=Inf)
	# ncd stands for number of chain dimensions
	# nsd stands for number of stored dimensions
	# this function returns a segmented complex with data for the first ncd
	# dimensions (the others are not included in the output); the output takes
	# the form of four length-ncd arrays

	nsd 	= 	length(dp)-1
	if ncd == Inf
		ncd  	= nsd
	end
	m 		= 	min(nsd,ncd)

	fvc     = Array{Array{Float64,1}}(undef,ncd)
	rvc     = Array{Array{Int64,1}}(undef,ncd)
	cpc     = Array{Array{Int64,1}}(undef,ncd)

	# remarks:
	# (a) #{cells of dimension ≤ (p   = k+1)} 				= dp[k+2]-1		= 	dp[p+1]-1
	# (b) #{cells of dimension ≤ (p-2 = (k-2)+1 = k-1)} 	= dp[k]-1 		=	dp[p-1]-1
	for p = 1:m
		rvc[p],cpc[p]   = 	copycolumnsubmatrix(rv,cp,cran(dp,p))
		rvc[p]			= 	rvc[p] .- (ec(dp,p-1,1) - 1)  # we subtract off the number of cells of dimension 2 less than those of interest, since starts the _faces_ of the cells of interest off at index 1; we get away with this in dimension 0 bc rvc[1] is the empty set
		fvc[p] 			=   convert(Array{Float64,1},fv[cran(dp,p)])
	end

	for p = (m+1):(ncd)
		rvc[p]		=	zeros(Int64,0)
		cpc[p]		=	ones(Int64,1)
		fvc[p] 		=	zeros(Int64,0)
	end

	if length(dp) 	>	ncd+1
		dpc			=	copy(dp[1:ncd+1])
	elseif length(dp)<  ncd+1
		dpc 		=	vcat(dp,fill(dp[end],ncd+1-length(dp))) # extend dp to the proper length
	else
		dpc 		=	copy(dp)
	end
	return rvc,cpc,fvc,dpc
end

function segmentedfilteredcomplex2unsegmentedfilteredcomplex(rv,cp,fv)
	numsd 					 	= 	length(fv)
	eulervec 					= 	zeros(Int64,numsd)
	for p 						= 	1:numsd
		eulervec[p]				= 	length(cp[p])-1
	end
	dp 							= 	eulervector2dimensionpattern(eulervec)

	cp 							= 	copy(cp)
	for 	p 					= 	2:numsd
		cp[p]					= 	(cp[p-1][end]-1).+cp[p]
		deleteat!(cp[p],1)
	end
	cp 							= 	cat(cp...,dims=1)

	rv 							= 	copy(rv)
	for 	p 					= 	2:numsd
		rv[p]					= 	rv[p] .+ (ec(dp,p-1,1) - 1)
	end
	rv 							= 	cat(rv...,dims=1)

	fv 							= 	cat(fv...,dims=1)
	return rv,cp,fv,dp
end


#=
version with 4 (integer array) non-keyword arguments
- 	rv:	row values
- 	cp: column pattern
- 	fv: filtration values
- 	dp: dimension pattern (satisfies cran(dp,j) = {cells of dimension j-1})

NB: the default value for maxdim has not been tested, and may cause errors;
the -3 accounts for (1) julia uses 1-indexed arrays, (2) to calculate
homology in dimension p, one must inspect the (p+1)-dimensional boundary
operator, (3) this operator should be givent the same treatment as those
that precede it ... generally this assumes that the next one up is at least
defined, even if it is trivial.
=#
function persistf2complex(	;
							rv 			= zeros(Int64,0),
							cp 			= zeros(Int64,0),
							fv 			= zeros(Int64,0),
							dp 			= zeros(Int64,0),
							dv 			= zeros(Int64,0),
							ev 			= zeros(Int64,0),
							maxdim		= [],
							numrad 		= Inf,
							maxrad 		= Inf,
							minrad 		= -Inf,
							record 		= "cyclerep",
							verbose		= false)

	if isempty(maxdim)
		if isempty(rv)
			ncd 			= 	1;
		elseif typeof(rv)  != 	Array{Int64,1}
			ncd 			= 	length(rv)-1 # we subtract one b/c it's usually convenient to have an extra empty array at the end of rv
		elseif !isempty(ev)
			ncd 			= 	length(ev)
		elseif !isempty(dv)
			ncd 			= 	dv[end]+1
			if !isempty(checkdv(rv,cp,dv))
				println("error: please check that the input operator is graded of degree 1")
				return
			# else
			# 	println("input operator is graded of degree 1")
			end
		elseif !isempty(dp)
			ncd 			= 	length(dp)-1
		else
			println("error: keyword input <rv> is a vector, and none of <dv>, <dp>, and <ev> is nonempty.")
		end
	else
		ncd 				= 	maxdim+2
	end


	if typeof(rv) == Array{Int64,1}
		if isempty(dp)
			if !isempty(ev)
				dp 				= 	eulervector2dimensionpattern(ev)
			end
			if !isempty(dv)
				dp 				= 	dimensionvalues2dimensionpattern(dv)
			end
		end
	end

	if !isempty(dp) # segment the complex if necessary
		rv,cp,fv,dp 	= 	unsegmentedfilteredcomplex2segmentedfilteredcomplex(rv,cp,fv,dp;ncd=ncd)
		if verbose
			if !pairwiseisequal([rv,cp,fv],under=length)
				println("Error: please check unsegmentedfilteredcomplex2segmentedfilteredcomplex")
			end
		end
	else # if it has not been defined by this point, define the dimension pattern
		ev 			= 	zeros(Int64,ncd)
		for p 			= 	1:ncd
			ev[p] 	= 	length(fv[p])
		end
		dp 				= 	eulervector2dimensionpattern(ev)
	end

	### Record the input parameters
	input = Dict(
		"model"			=> "complex",
		# "version" 		=> Pkg.installed("Eirene"),
		"date"			=> string(Dates.Date(now())),
		"time"			=> string(Dates.Time(now())),
		"maxdim"		=> ncd-2,
		"record" 		=> record
		)

	if typeof(fv) == Array{Float64,1}
		# println("nonerror message: TYPE IS FLOATING ARRAY!")
		ocg,ocg2rad			=   trueordercanonicalform(
								fv,
								firstval=1,
								factor=true,
								rev=true)
	else
	 	ocg,ocg2rad			=   trueordercanonicalform(
								cat(fv...,dims=1),
								firstval=1,
								factor=true,
								rev=true)
	end

	ocg 					= 	segmentarray(ocg,dp)

	### Perform the persistence computation
	trv,tcp,plo,phi,tid =
	persistf2_core_cell(
		rv,
		cp,
		ocg;
		maxsd = ncd,
		record=record,
		verbose=false,
		prepairs = convert(Array{Array{Int64,1},1},fill(Array{Int64}(undef,0),ncd))
		)

	### Create the dictionary that stores all relevant data
	D = Dict{String,Any}(
		"rv" 		=> rv,
		"cp" 		=> cp,
		"grain"		=> ocg,
		"trv" 		=> trv,
		"tcp" 		=> tcp,
		"tid" 		=> tid,
		"plo" 		=> plo,
		"phi" 		=> phi,
		"ocg2rad" 	=> ocg2rad,
		"input"		=> input
		)

	#### Store generators
	#gc()
	if record == "all" || record == "cyclerep"
 		unpack!(D)
	end
	#gc()
	if record == "cyclerep"
		delete!(D,"trv")
		delete!(D,"tcp")
		delete!(D,"Lrv")
		delete!(D,"Lcp")
		delete!(D,"Rrv")
		delete!(D,"Rcp")
		delete!(D,"Lirv")
		delete!(D,"Licp")
		delete!(D,"prepairs")
	end

	return D
end


function parse_perseusdistmat(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,:])
	minrad = convert(Float64,A[2,1])
	maxrad = minrad+convert(Float64,A[2,2])*convert(Float64,A[2,3])
	maxsd = convert(Int64,A[2,4])+1
	return s,minrad,maxrad,maxsd
end

function parse_perseusbrips(filename)
	A = readdlm(filename)
	s = convert(Array{Float64},A[3:end,1:end-1])
	minrad = 0
	maxrad = convert(Array{Float64},A[2,1:3])
	maxrad = prod(maxrad)
	return s,minrad,maxrad
end

function boundarylimit_Lionly(trv,tcp,tid,sd;verbose=false)
	if isempty(tid[sd])
		Lirv = zeros(Int64,0)
		Licp = ones(Int64,1)
	else
		tsize			= length(tcp[sd])-1
		Lrv 			= copy(trv[sd])
		lowtranslator 	= zeros(Int64,maximum(tid[sd]))
		lowtranslator[tid[sd]] = 1:length(tid[sd])
		yafterx!(lowtranslator,Lrv)
		Lirv,Licp   	= morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv,tcp[sd])
		Lirv,Licp   	= transposeLighter(Lirv,Licp,tsize)
	end
	return Lirv,Licp
end

#=
-	L, Li, and R are all indexed by tid, just like (trv,tcp)
- 	L (not Li) is the transpose of (trv,tcp)
- 	up to permutation and deletion of some zero rows, RLM is the fundamental
 	cycle matrix of the d-dimensional boundary operator with respect
	to the doubly-minimal basis we have constructed, where M the submatrix of
	the boundary with rows indexed by tid.
=#
function boundarylimit_core(brv,bcp,trv,tcp,tid,numl,nump,numnlpl)
	lowtranslator = zeros(Int64,numl)
	lowtranslator[tid] = 1:numnlpl
	trv = copy(trv)
	yafterx!(lowtranslator,trv)
	#################################################################################
	# if !issorted(tcp)
	# 	println("error: tcp is not sorted")
	# 	sleep(2)
	# end
	# if length(trv) != tcp[end]-1
	# 	println("error: tcp appears to end in the wrong place")
	# 	sleep(2)
	# end
	# if tcp[1] != 1
	# 	println("error: tcp does not begin with 1")
	# 	sleep(2)
	# end
	#################################################################################
	Lirv,Licp	= morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(trv,tcp)
	Lirv,Licp   = transposeLighter(Lirv,Licp,numnlpl)
	Lrv,Lcp 	= transposeLighter(trv,tcp,numnlpl)
	brv,bcp 	= spmmF2silentLeft(Lrv,Lcp,brv,bcp,numnlpl)
	Rrv,Rcp 	= morseInverseF2orderedColsUnsortedRowsInSilentOut(brv,bcp)
	return Lrv,Lcp,Lirv,Licp,Rrv,Rcp
end

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
	numl 		= length(farfaces[sd-1])
	nump 		= length(phi[sd])
	numnlpl 	= numl-length(plo[sd-1])
	brv,bcp 	= maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
	return 	boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
	numl 		= length(cp[sd-1])-1
	nump 		= length(phi[sd])
	numnlpl 	= numl-length(phi[sd-1])
	brv,bcp		= maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
	Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
	return 	boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit(D::Dict,sd)
	# special case sd = 1 may possibly be degenerate
	trv = D["trv"];tcp=D["tcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
	if haskey(D,"farfaces")
		farfaces = D["farfaces"];firstv = D["firstv"]
		return boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
	else
		rv = D["rv"];cp = D["cp"]
		# Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
		return boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
	end
end

function maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
	numl = length(farfaces[sd-1])
	nump = length(phi[sd])
	numnlpl = numl-length(plo[sd-1])

	lowtranslator = zeros(Int64,numl)

	lowtranslator[tid[sd]] = 1:numnlpl
	brv = ff2aflight(farfaces,firstv,sd,phi[sd])
	brv = reshape(brv,length(brv))
	bcp = convert(Array{Int64,1},1:sd:(nump*sd+1))
	supportedmatrix!(brv,bcp,tid[sd],1:nump,numl)
	yafterx!(lowtranslator,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. too long
	return brv,bcp
end

function maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
	numl = length(cp[sd-1])-1
	nump = length(phi[sd])
	numnlpl = numl-length(plo[sd-1])

	if isempty(phi[sd])
		brv 	= 	zeros(Int64,0)
		bcp 	=   ones(Int64,numnlpl+1)
		return 		brv,bcp
	end

	lowtranslator = zeros(Int64,numl)
	lowtranslator[tid[sd]] = 1:numnlpl
	dummy0 = zeros(Int64,0)

	brv,dummy1,bcp,dummy2 = stackedsubmatrices(rv[sd],cp[sd],tid[sd],dummy0,phi[sd],numl)
	yafterx!(lowtranslator,brv)
	append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. longer than the # of columns
	return brv,bcp
end

function unpack!(D::Dict)
	# l = length(D["grain"])
	maxsd = D["input"]["maxdim"]+2 # grain2maxsd(D["grain"])

	Lirv = Array{Array{Int64,1},1}(undef,maxsd);  Licp = Array{Array{Int64,1},1}(undef,maxsd)
	Lrv  = Array{Array{Int64,1},1}(undef,maxsd);  Lcp  = Array{Array{Int64,1},1}(undef,maxsd)
	Rrv  = Array{Array{Int64,1},1}(undef,maxsd);  Rcp  = Array{Array{Int64,1},1}(undef,maxsd)

	if 	D["input"]["record"] == "all"
		N 	= 	maxsd
	else
		N 	= 	maxsd-1
		if maxsd >= 2
			Lirv[maxsd],Licp[maxsd] = boundarylimit_Lionly(D["trv"],D["tcp"],D["tid"],maxsd)
		elseif maxsd == 1
			Lirv[maxsd] = Array{Int64,1}(undef,0)
			Licp[maxsd] = ones(Int64,1)
		end
		Lrv[maxsd]=Array{Int64,1}(undef,0)
		Lcp[maxsd]=Array{Int64,1}(undef,0)
		Rrv[maxsd]=Array{Int64,1}(undef,0)
		Rcp[maxsd]=Array{Int64,1}(undef,0)
	end

	for i = 2:N
		Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = boundarylimit(D,i)
		if isempty(Lcp[i])
			println()
			println("ERROR MESSAGE IN unpack!: Lcp[i] = 0 and i = $(i)")
		end
	end

	Lirv[1]		=Array{Int64}(undef,0)
	Lrv[1]		=Array{Int64}(undef,0);
	Rrv[1]		=Array{Int64}(undef,0)
	Licp[1]		=ones(Int64,1)
	Lcp[1]		=ones(Int64,1)
	Rcp[1]		=ones(Int64,1)

	D["Lrv"] = Lrv
	D["Lcp"] = Lcp
	D["Lirv"]= Lirv
	D["Licp"]= Licp
	D["Rrv"] = Rrv
	D["Rcp"] = Rcp

	D["cyclerep"] = fill(Array{Array{Int64,1},1}(undef,0),maxsd)  # for each dimension one stores an array of arrays

	for i = 2:maxsd
		dim = i-2
		m = nnzbars(D,dim=dim)
		cyclenames = barname2cyclename(D,1:m,dim=dim)
		D["cyclerep"][i] = getcycle(D,cyclenames,dim=dim)
	end
	return
end

##########################################################################################

####	INVERSION

##########################################################################################

function morseInverseF2orderedColsUnsortedRowsInSilentOut(Arowval::Array{Tv,1},Acolptr::Array{Tv,1}) where Tv<:Integer
	mA = length(Acolptr)-1
	colptrA = Acolptr
	rowvalA = Arowval
    preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(undef,mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}(undef,mA)
	rowvalCj = Array{Bool}(undef,mA)
	rowvalC = Array{Tv}(undef,mA)
    totalrowscounter = 0
    onepast = 0
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array{Int64}(undef,preallocationIncrement))
		end
		if colptrA[i]+1 == colptrA[i+1]
			colptrC[i+1] = colptrC[i]
		elseif colptrA[i]+2==colptrA[i+1]
			if Arowval[colptrA[i]]<i
				k = Arowval[colptrA[i]]
			else
				k = Arowval[colptrA[i]+1]
			end
			ccap = colptrC[i]+colptrC[k+1]-colptrC[k]-1
			rowvalC[colptrC[i]:ccap]= rowvalC[colptrC[k]:(colptrC[k+1]-1)]
			rowvalC[ccap+1]=k
			colptrC[i+1]=ccap+2
		else
			eyerange = cran(Acolptr,i)
			newRowsCounter = length(eyerange)
			for l=1:newRowsCounter
				ll = Arowval[eyerange[l]]
				rowList[l] = ll
				rowSupp[ll] = i
				rowvalCj[ll]=true
			end
			rowvalCj[i] = false ## note we have to make this correction
			for jp in eyerange
				j = rowvalA[jp]
				if j < i
					for kp in colptrC[j]:(colptrC[j+1] - 1)
						k = rowvalC[kp]
						if rowSupp[k] != i
							rowSupp[k] = i
							newRowsCounter +=1
							rowList[newRowsCounter] = k
							rowvalCj[k] = true
						else
							rowvalCj[k] = !rowvalCj[k]
						end
					end
				end
			end
			marker = colptrC[i]
			for l = 1:newRowsCounter
				if rowvalCj[rowList[l]]
					rowvalC[marker]=rowList[l]
					marker+=1
				end
			end
			colptrC[i+1]=marker
		end

	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Arowval::Array{Tv,1},Acolptr::Array{Tv,1}) where Tv<:Integer
	mA = length(Acolptr)-1
	colptrA = Acolptr
	rowvalA = Arowval
    preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(undef,mA+1); colptrC[1]=1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}(undef,mA)
	rowvalCj = Array{Bool}(undef,mA)
	rowvalC = Array{Tv}(undef,mA)
    totalrowscounter = 0
    onepast = 0
	for i in 1:mA
		if colptrC[i]+mA > length(rowvalC)+1
			append!(rowvalC, Array{Int64}(undef,preallocationIncrement))
		end
		if colptrA[i] == colptrA[i+1]
			colptrC[i+1] = colptrC[i]
		elseif colptrA[i]+1==colptrA[i+1]
			k = Arowval[colptrA[i]]
			ccap = colptrC[i]+colptrC[k+1]-colptrC[k]
			rowvalC[colptrC[i]:(ccap-1)]= crows(colptrC,rowvalC,k)
			rowvalC[ccap]=k
			colptrC[i+1]=ccap+1
		else
			eyerange = cran(Acolptr,i)
			newRowsCounter = length(eyerange)
			for l=1:newRowsCounter
				ll = Arowval[eyerange[l]]
				rowList[l] = ll
				rowSupp[ll] = i
				rowvalCj[ll]=true
			end
			for jp in eyerange
				for kp in cran(colptrC,rowvalA[jp])
					k = rowvalC[kp]
					if rowSupp[k] != i
						rowSupp[k] = i
						newRowsCounter +=1
						rowList[newRowsCounter] = k
						rowvalCj[k] = true
					else
						rowvalCj[k] = !rowvalCj[k]
					end
				end
			end
			marker = colptrC[i]
			for l = 1:newRowsCounter
				if rowvalCj[rowList[l]]
					rowvalC[marker]=rowList[l]
					marker+=1
				end
			end
			colptrC[i+1]=marker
		end
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

##########################################################################################

####	MULTIPLICATION

##########################################################################################

function spmmF2(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am) where Tv<:Integer
    mA = Am
    nB = length(Bcolptr)-1
    rowvalA = Arowval; colptrA = Acolptr
    rowvalB = Browval; colptrB = Bcolptr
    preallocationIncrement = colptrA[end]+colptrB[end]

	colptrC = Array{Tv}(undef,nB+1)
    colptrC[1] = 1
	rowSupp = zeros(Tv, mA)
	rowList = Array{Tv}(undef,mA)
	rowvalCj = Array{Bool}(undef,mA)
	rowvalC = Array{Tv}(undef,preallocationIncrement)
	for i in 1:nB
		newrowscounter = 0
		for jp in colptrB[i]:(colptrB[i+1] - 1)
			j = rowvalB[jp]
			for kp in colptrA[j]:(colptrA[j+1] - 1)
				k = rowvalA[kp]
				if rowSupp[k] != i
					rowSupp[k] = i
					newrowscounter +=1
					rowList[newrowscounter] = k
					rowvalCj[k] = true
				else
					rowvalCj[k] = !rowvalCj[k]
				end
			end
		end
		nzRows = findall(rowvalCj[rowList[1:newrowscounter]])
		colptrC[i+1] = colptrC[i]+length(nzRows)

		if colptrC[i+1] > length(rowvalC)+1
			append!(rowvalC, Array{Int}(undef,preallocationIncrement))
		end
		rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function spmmF2silentLeft(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am) where Tv<:Integer
    mA = Am
    nB = length(Bcolptr)-1
    rowvalA = Arowval; colptrA = Acolptr
    rowvalB = Browval; colptrB = Bcolptr
    preallocationIncrement = colptrA[end]+colptrB[end]

	colptrC = Array{Tv}(undef,nB+1)
    colptrC[1] = 1
	rowSupp = zeros(Tv,mA)
	rowList = Array{Tv}(undef,mA)
	rowvalCj = Array{Bool}(undef,mA)
	rowvalC = Array{Tv}(undef,preallocationIncrement)
	for i in 1:nB
		newrowscounter = 0
		eyerange = cran(Bcolptr,i)
		newrowscounter = length(eyerange)
		for l=1:newrowscounter
			ll = rowvalB[eyerange[l]]
			rowList[l] = ll
			rowSupp[ll] = i
			rowvalCj[ll]=true
		end
		for jp in colptrB[i]:(colptrB[i+1] - 1)
			j = rowvalB[jp]
			for kp in colptrA[j]:(colptrA[j+1] - 1)
				k = rowvalA[kp]
				if rowSupp[k] != i
					rowSupp[k] = i
					newrowscounter +=1
					rowList[newrowscounter] = k
					rowvalCj[k] = true
				else
					rowvalCj[k] = !rowvalCj[k]
				end
			end
		end
		nzRows = findall(rowvalCj[rowList[1:newrowscounter]])
		colptrC[i+1] = colptrC[i]+length(nzRows)

		if colptrC[i+1] > length(rowvalC)+1
			append!(rowvalC, Array{Int}(undef,preallocationIncrement))
		end
		rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
	end
	deleteat!(rowvalC,colptrC[end]:length(rowvalC))
	return rowvalC, colptrC
end

function blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement)
	extend!(Mcp,Dn+1)
	deleteat!(Mcp,(Dn+2):length(Mcp))
	Mcp[1]=1

	rowSupp = zeros(Int64, Dm)
	rowList = Array{Int64}(undef, Dm)
	rowvalCj = BitArray(undef,Dm)
	for i in 1:Dn
		if length(Mrv)<Mcp[i]+Dm
			extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
		end
		eyerange = cran(Dcp,i)
		newrowscounter = 0
		for kp = eyerange
			row = Drv[kp]
			rowSupp[row]=i
			newrowscounter+=1
			rowList[newrowscounter]=row
			rowvalCj[row] = true
		end
		if Bcp[i]<Bcp[Dn+1]
			for jp in Bcp[i]:(Bcp[i+1] - 1)
				j = Brv[jp]
				for kp in Ccp[j]:(Ccp[j+1] - 1)
					row = Crv[kp]
					if rowSupp[row] != i
						rowSupp[row] = i
						newrowscounter +=1
						rowList[newrowscounter] = row
						rowvalCj[row] = true
					else
						rowvalCj[row] = !rowvalCj[row]
					end
				end
			end
		end
		ii = i+1
		Mcp[ii]=Mcp[i]
		for k = 1:newrowscounter
			row = rowList[k]
			if rowvalCj[row]
				Mrv[Mcp[ii]]=row
				Mcp[ii]+=1
			end
		end
	end
end

function blockprodsumsilenticolsleftWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement,col2silenti)
	extend!(Mcp,Dn+1)
	deleteat!(Mcp,(Dn+2):length(Mcp))
	Mcp[1]=1

	rowSupp = zeros(Int64, Dm)
	rowList = Array{Int64}(undef,Dm)
	rowvalCj = BitArray(undef,Dm)
	for i in 1:Dn
		if length(Mrv)<Mcp[i]+Dm
			extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
		end
		eyerange = cran(Dcp,i)
		newrowscounter = 0
		for kp = eyerange
			row = Drv[kp]
			rowSupp[row]=i
			newrowscounter+=1
			rowList[newrowscounter]=row
			rowvalCj[row] = true
		end
		for jp in Bcp[i]:(Bcp[i+1] - 1)
			j = Brv[jp]
			row = col2silenti[j]
			if rowSupp[row] != i
				rowSupp[row] = i
				newrowscounter +=1
				rowList[newrowscounter] = row
				rowvalCj[row] = true
			else
				rowvalCj[row] = !rowvalCj[row]
			end
		end
		for jp in Bcp[i]:(Bcp[i+1] - 1)
			j = Brv[jp]
			for kp in Ccp[j]:(Ccp[j+1] - 1)
				row = Crv[kp]
				if rowSupp[row] != i
					rowSupp[row] = i
					newrowscounter +=1
					rowList[newrowscounter] = row
					rowvalCj[row] = true
				else
					rowvalCj[row] = !rowvalCj[row]
				end
			end
		end
		ii = i+1
		Mcp[ii]=Mcp[i]
		for k = 1:newrowscounter
			row = rowList[k]
			if rowvalCj[row]
				Mrv[Mcp[ii]]=row
				Mcp[ii]+=1
			end
		end
	end
end

##########################################################################################

####	COPY, SHIFT, INDEX, AND SLICE

##########################################################################################

function cran(A::SparseMatrixCSC,j)
	return A.colptr[j]:(A.colptr[j+1]-1)
end

function cran(colptr::Array,j::Int64)
	return colptr[j]:(colptr[j+1]-1)
end

#	added 12/27/2017
# 	may be tested via function <testcran>
function cran(colptr::Array,J::Array{Int64,1})
	m = nval(colptr,J)
	v = zeros(Int64,m)
	c = 0
	for p=1:length(J)
		k = nval(colptr,J[p])
		v[c+1:c+k]=cran(colptr,J[p])
		c += k
	end
	return v
end

#	added 12/27/2017
# 	may be tested via function <testcran>
function cran(colptr::Array,J::UnitRange{Int64})
	m = nval(colptr,J)
	v = zeros(Int64,m)
	c = 0
	for p=1:length(J)
		k = nval(colptr,J[p])
		v[c+1:c+k]=cran(colptr,J[p])
		c += k
	end
	return v
end


function cran(colptr::UnitRange,j)
	return colptr[j]
end

#	added 12/27/2017
"""
	nval(colptr,j::Int64)

For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column `j`.
"""
function nval(colptr,j::Int64)
	return colptr[j+1]-colptr[j]
end

#	added 12/27/2017
"""
	nval(colptr,J)

For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column j, for each j in `J`.
"""
function nval(colptr,J)
	c = 0
	for p = 1:length(J)
		c += nval(colptr,J[p])
	end
	return c
end

function crows(A::SparseMatrixCSC,j)
	return A.rowval[cran(A,j)]
end

function crows(colptr::Array,rowval::Array,j)
	return rowval[cran(colptr,j)]
end

function extend!(x::Array{Tv,1},n::Integer) where Tv
	if length(x)<n
		append!(x,Array{Tv}(undef,n-length(x)))
	end
end

function copycolumnsubmatrix(Arv::Array{Tv,1},Acp,columnindices) where Tv<:Integer
	allocationspace = 0
	for j in columnindices
		allocationspace+= Acp[j+1]-Acp[j]
	end
	Brv = Array{Tv}(undef,allocationspace)
	Bcp = Array{Int64}(undef,length(columnindices)+1)
	Bcp[1]=1
	for jp = 1:length(columnindices)
		j = columnindices[jp]
		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
	end
	return Brv,Bcp
end

# NB: the following was a facsimile of the preceding funciton, written at a time
# when it was regarded as a good style in Julia to be extremely specific with
# regard to input type.
# function copycolumnsubmatrix{Tv<:Integer}(Arv::Array{Tv,1},Acp::Array{Tv,1},columnindices::UnitRange{Int64})
# 	allocationspace = 0
# 	for j in columnindices
# 		allocationspace+= Acp[j+1]-Acp[j]
# 	end
# 	Brv = Array{Tv}(undef,allocationspace)
# 	Bcp = Array{Tv}(undef,length(columnindices)+1)
# 	Bcp[1]=1
# 	for jp = 1:length(columnindices)
# 		j = columnindices[jp]
# 		Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
# 		Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
# 	end
# 	return Brv,Bcp
# end

function copycolind2colind!(rowvalA::Array{Tv,1},colptrA::Array{Tv,1},columnindices,rowvalB::Array{Tv,1},colptrB::Array{Tv,1},startingDestination::Integer,growthIncrement::Integer)  where Tv<:Integer
	numnewrows = 0
	for col in columnindices
		numnewrows+=colptrA[col+1]-colptrA[col]
	end
	if length(rowvalB)<colptrB[startingDestination]+numnewrows-1
		append!(rowvalB,Array{Tv}(undef,numnewrows+growthIncrement))
	end
	if length(colptrB)<startingDestination+length(columnindices)
		append!(colptrB,Array{Tv}(undef,startingDestination+length(columnindices)))
	end
	colptrB[1]=1
	for i = 1:length(columnindices)
		k = startingDestination+i #the index of colptr pointing to end of this new column
		col = columnindices[i]
		colptrB[k]=colptrB[k-1]+colptrA[col+1]-colptrA[col]
		rowvalB[colptrB[k-1]:(colptrB[k]-1)]=rowvalA[colptrA[col]:(colptrA[col+1]-1)]
	end
end


# colsinorder must be in sorted order
function supportedmatrix!(Mrowval::Array{Tv},Mcolptr::Array{Tv,1},rows1,colsinorder,Mm::Tv) where Tv<:Integer
	n = length(colsinorder)
	suppcol1 = falses(Mm)
	suppcol1[rows1].=true
	cpHolder = 1
	nz1 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,colsinorder[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
				Mrowval[nz1]=i
			end
		end
		Mcolptr[jp]=cpHolder
		cpHolder = nz1+1
	end
	Mcolptr[n+1] = cpHolder
	deleteat!(Mcolptr,(n+2):length(Mcolptr))
	deleteat!(Mrowval,Mcolptr[end]:length(Mrowval))
end

function stackedsubmatrices(
	Mrowval,#::Array{Tv,1},
	Mcolptr,#::Array{Tv,1},
	rows1,#::Array{Tv,1},
	rows2,#::Array{Tv,1},
	cols,#::Array{Tv,1},
	Mm::Tv)  where Tv<:Integer

	n = length(cols)
	suppcol1 = falses(Mm)
	suppcol2 = falses(Mm)
	suppcol1[rows1].=true
	suppcol2[rows2].=true
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
			elseif suppcol2[i]
				nz2+=1
			end
		end
	end
	rv1 = Array{Tv}(undef,nz1)
	rv2 = Array{Tv}(undef,nz2)
	cp1 = Array{Tv}(undef,n+1); cp1[1]=1
	cp2 = Array{Tv}(undef,n+1); cp2[1]=1
	nz1 = 0; nz2 = 0
	for jp = 1:n
		for ip in cran(Mcolptr,cols[jp])
			i = Mrowval[ip]
			if suppcol1[i]
				nz1+=1
				rv1[nz1]=i
			elseif suppcol2[i]
				nz2+=1
				rv2[nz2]=i
			end
		end
		cp1[jp+1] = nz1+1
		cp2[jp+1] = nz2+1
	end
	return rv1,rv2,cp1,cp2
end

##########################################################################################

####	TRANSPOSITION OPERATORS

##########################################################################################

function transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Am) where Ti
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    # Attach destination matrix
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(undef,Am+1)
    Crowval = Array{Ti}(undef,Annz)
    # Compute the column counts of C and store them shifted forward by one in
	# Ccolptr
    Ccolptr[1:end] .= 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and
	# Cnzval, tracking write positions in Ccolptr
    for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr
	# shift, but the first colptr remains incorrect
    Ccolptr[1] = 1

	return Crowval, Ccolptr
end

function transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Anzval::Array{Tv},Am::Integer) where {Tv, Ti}
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(undef,Am+1)
    Crowval = Array{Ti}(undef,Annz)
    Cnzval = Array{Tv}(undef,Annz)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] .= 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and Cnzval,
    # tracking write positions in Ccolptr
    @inbounds for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Cnzval[Ck] = Anzval[Ak]
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1
	return Crowval, Ccolptr, Cnzval
end

#=
-	Returns the sparse equivalent of A[rows,cols]'.
- 	The rows of the output array will be listed in the order that they appear in
	input vector <cols>.
- 	Does not assume that the entries of <rows> appear in sorted order.
- 	I *believe* that repeated rows and columns are allowed.
=#
function transposeLighter_submatrix(Arowval::Array{Ti},Acolptr::Array{Ti},Am;rows = 1:Am,cols = 1:length(Acolptr)-1) where Ti
	if rows == 1:Am && cols == 1:(length(Acolptr)-1)
		Crowval, Ccolptr = transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Am)
		return Crowval, Ccolptr
	end
    # Attach destination matrix
    Cm = length(cols)
    Cn = length(rows)
    Ccolptr = Array{Ti}(undef,Cn+1)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] .= 0
	rs = rowsupportsum(Arowval,Acolptr,Am,cols)
	for i = 1:Cn
	    Ccolptr[i+1] = rs[rows[i]]
	end
	Cnnz = sum(Ccolptr)
    Crowval = Array{Ti}(undef,Cnnz)
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and Cnzval,
    # tracking write positions in Ccolptr
    rowtranslator = zeros(Int64,Am)
    rowtranslator[rows] = 1:Cn
    @inbounds for Cj in cols
        for Ak in cran(Acolptr,Cj)
            Ai = rowtranslator[Arowval[Ak]]
            if Ai > 0
	            Ck = Ccolptr[Ai+1]
	            Crowval[Ck] = Cj
    	        Ccolptr[Ai+1] += 1
    	    end
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1
	return Crowval, Ccolptr
end

#=
Accepts an mxn integer array, M.  To M we implicitly associate an array N,
as follows. If M has no nonzero entries, then N is the zeron-one array such that
supp(N[:,i]) = M[:,i].  Here we regard M[:,i] as a set.  I *believe* one ignores
duplicate entries, but have not checked. A similar interpretation holds when M
has zero entries - one simply discards the zeros.
Let S denote the support of N, r denote row02row1translator, and c denote
col02col1translator.  This function returns the data specifying a sparse
zero-one matrix whose support is
{ [c[j],r[i]] : [i,j] \in N and neither c[j] nor r[i] are zero }.
=#
function presparsefull2unsortedsparsetranspose(
	M::Array{Tv,2},
	row02row1translator,
	col02col1translator;
	verbose::Bool=false) where Tv<:Integer

	Mm,Mn = size(M)

	if Mn == 0
		rowval1 = Array{Int64}(undef,0)
		if isempty(row02row1translator)
			colptr1 = ones(Int64,1)
		else
			colptr1 = ones(Int64,1+maximum(row02row1translator))
		end
		return rowval1,colptr1,Mn
	end

	for i = 1:(Mm*Mn)
		if M[i]>0
			M[i] = row02row1translator[M[i]]
		end
	end
	m0 = maximum(row02row1translator)
	rowcounter = zeros(Int64,m0)
	for k in M  #allow for zero values
		if k>0
			rowcounter[k]+=1
		end
	end
	colptr1 = Array{Int64}(undef,m0+1)
	colptr1[1]=1
	for i = 1:m0
		colptr1[i+1]=colptr1[i]+rowcounter[i]
	end
	rowval1 = Array{Int64}(undef,colptr1[end]-1)
	placer = copy(colptr1)
	if verbose
		coverageTestVec = trues(colptr1[end]-1)
	end
	for j = 1:Mn
		jval = col02col1translator[j]
		for i = 1:Mm
			row = M[i,j]
			if row > 0
				if verbose
					coverageTestVec[placer[row]]=false
				end
				rowval1[placer[row]]=jval
				placer[row]+=1
			end
		end
	end
	if verbose
		if any(coverageTestVec)
			print("please refer to the coverageTestVec")
			sleep(4)
		end
		println([length(rowval1) "length(output of presparsefull2unsortedsparsetranspose)"])
	end
	#gc()
	return rowval1,colptr1,Mn
end

##########################################################################################

####	SHAPE GENERATORS

##########################################################################################

function noisycircle()
	theta = 1:100
	theta = theta*2*pi/100
	x = cos.(theta)
	y = sin.(theta)
	pcloud = hcat(x,y)'
	pcloud = hcat(pcloud,pcloud,pcloud)
	pcloud = pcloud + 0.3*rand(2,300)
	return pcloud
end

function noisycircle3()
	theta = 1:100
	theta = theta*2*pi/100
	x = cos.(theta)
	y = sin.(theta)
	z = zeros(100)
	pcloud = hcat(x,y,z)'
	pcloud = hcat(pcloud,pcloud,pcloud)
	pcloud = pcloud + 0.3*rand(3,300)
	return pcloud
end

function torus(;m = 100,n = 50,mrad=1,nrad = 0.2)
	theta = (1:m)*2*pi/m;
	torus1 = mrad*repeat(cos.(theta),1,n)
	torus2 = mrad*repeat(sin.(theta),1,n)
	torus3 = nrad*repeat(sin.(2*pi*(1:n)./n),1,m)'
	torus4 = repeat(cos.(2*pi*(1:n)./n),1,m)'
	torus1 = torus1+nrad*torus1.*torus4;
	torus2 = torus2+nrad*torus2.*torus4;
	return hcat(torus1[:],torus2[:],torus3[:])'
end

function noisytorus(;m = 100,n=50,mrad=1,nrad = 0.2,noiserad= 0.5*nrad)
	pcloud = torus(m=m, n=n,nrad = nrad,mrad=mrad)
	pcloud = pcloud + noiserad*rand(size(pcloud,1),size(pcloud,2))
end

function sphere()
	halfnumrings = 40
	latitd = pi*(-halfnumrings:halfnumrings)/(2*halfnumrings)

	x = Array{Float64}(undef,0)
	y = Array{Float64}(undef,0)
	z = Array{Float64}(undef,0)

	for t in latitd
		numpts = round.(Int64,40*cos.(t))
		theta = 2*pi*(1:numpts)/numpts
		append!(x,cos.(t)*cos.(theta))
		append!(y,cos.(t)*sin.(theta))
		append!(z,fill(sin.(t),numpts))
	end
	return x,y,z
end

function matchingcomplex_symmat(m,n)
	#### 1-A-I, where A is the adjacency matrix of the matching complex M(m,n)
	#### and I is identity
	combs = collect(combinations(Array(1:n),m))  #an array of arrays
	numcombs = binom(n,m)
	numpointedcombs = binom(n-1,m-1)
	symmat = zeros(Int64,numcombs,numcombs)
	indexvec = zeros(Int64,numpointedcombs)
	marker = 0

	for i = 1:n
		marker = 0
		for j = 1:numcombs
			if i in combs[j]
				marker+=1
				indexvec[marker]=j
			end
		end
		if marker!=length(indexvec) || marker!=numpointedcombs
			print("please check construction algorithm")
			sleep(10)
		end
		symmat[indexvec,indexvec] = 1
	end
	for i = 1:numcombs
		symmat[i,i]=0
	end
	return symmat
end

function chessboardcomplex_symmat(;numrows=3,numcols=4)
	#### 1-A-I, where A is the adjacency matrix of the chessboard complex C(m,n)
	#### and I is identity
	m 		= numrows
	n 		= numcols
	numrooks = m*n
	symmat = ones(Int64,numrooks,numrooks)
	for i = 1:m
		for j = 1:n
			rooknum1 = (i-1)*m+j
			for ii = 1:m
				for jj = 1:n
					if ii != i && jj !=j
						rooknum2 = (ii-1)*m+jj
						symmat[rooknum1,rooknum2] = 0
						symmat[rooknum2,rooknum1] = 0
					end
				end
			end
		end
	end
	for k = 1:numrooks
		symmat[k,k] = 0
	end
	return symmat
end

function plane2torus(A)
	theta1 = A[:,1]
	theta2 = A[:,2]
	x = cos.(theta1)
	y = sin.(theta1)
	z = 0.25*sin.(theta2)
	alpha = 1-0.25*cos.(theta2)
	x = x[:]'
	y = y[:]'
	z = z[:]'
	alpha = alpha[:]'
	x = alpha.*x
	y = alpha.*y
	return vcat(x,y,z)
end

function latlon2euc(A;model = "pc")
	if model == "pc"
		theta1 = A[1,:]
		theta2 = A[2,:]
	elseif model == "points"
		theta1 = A[:,1]'
		theta2 = A[:,2]'
	end
	x = cosd.(theta2)
	y = sind.(theta2)
	z = sind.(theta1)
	w = cosd.(theta1)
	x = w.*x
	y = w.*y
	return vcat(x,y,z)
end


function zerodrandmat(n)
	# input: 	an integer n
	# output: 	symmetric matrix with zeros on the diagonal and iid uniform entries elsewhere
	x = rand(n,n)
	for p = 1:n
		x[p,p]=0
	end
	x = (x+x')/2
	return x
end

function eirenefilepath(filedescription)
	if 		filedescription 	== 	"simplecity"
			return joinpath(@__DIR__,"examples/simplemapscitydata.csv")
	elseif 	filedescription 	== 	"noisycircle"
			return joinpath(@__DIR__,"examples/noisycircle.csv")
	elseif 	filedescription 	== 	"noisytorus"
			return joinpath(@__DIR__,"examples/noisytorus.csv")
	end
end

##########################################################################################

####	MATRIX WEIGHTS AND FORMATTING

##########################################################################################

#=
The output of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil(N;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)
	S 	= 	copy(N) # NB: this is probably an important step; several pernicious problems in the development of rounding procedures turned out to be rooted in accidental rewriting of matrix entries
	return minmaxceil!(S;minrad=minrad,maxrad=maxrad,numrad=numrad)
end

#=
NB THIS FUNCTION MODIFIES ITS INPUT

The result of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil!(S;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)

	if 	(minrad == Inf) || (maxrad == -Inf)
		return fill(Inf,size(S)...)
	end

	if minrad == "minedge"
		minrad = minimum(offdiagmin(S))
	end

	S[S.<minrad] 	.= 	minrad
	S[S.>maxrad] 	.= 	Inf

	fi 	= 	(LinearIndices(S))[findall(isfinite,S)] # stands for finite indices
	fv 	= 	S[fi]

	if isempty(fv)
		return S
	end

	if minrad 	== 	-Inf
		minrad 	=  	minimum(fv)
	end
	if maxrad 	== 	Inf
		maxrad 	= 	maximum(fv)
	end
	if numrad 	== 	1
		S[fi]	.= 	maxrad
		return 	S
	end
	if numrad 	== 	Inf
		return S
	end

	ran 		= 	range(minrad,maxrad,length=numrad)
	ran 		= 	Array{Float64}(ran)
	ran[end] 	= 	maxrad

	fvr 		= 	ceil2grid(fv,ran)

	S[fi]		= 	fvr

	return 		S
end

# ran should be an array in sorted order, with ran[end] >= maximum(A)
function ceil2grid(A,ran)
	if 	ran 			== 		"all"
		return A
	end
	B 					= 		copy(A)
	for j 				= 		1:length(A)
		post 			= 		1
		while ran[post] < A[j]
			post+=1
		end
		B[j] 			= 		ran[post]
	end
	return B
end


#=
minimum(deleteat(S[:,i],i))
=#
function offdiagmin(S,i)
	if i == 1
		return(minimum(S[2:end,i]))
	elseif i == size(S,1)
		return minimum(S[1:end-1,i])
	else
		return min(minimum(S[1:i-1,i]),minimum(S[i+1:end,i]))
	end
end


function offdiagmean(S;defaultvalue=[])
	m,n 	= 	size(S)
	if m 	!= 	n
		println()
		println("error in <offdiagmean>: input matrix should be square")
	end
	if isempty(defaultvalue)
		println()
		println("warning: no defaulvalue was set for <offdiagmin>; this parameter has been set to 0")
		defaultvalue 	= 	0
	end
	if m 		== 	1
		return 	defaultvalue
	end
	mu 			= 	zeros(m)
	for j 		= 	1:m
		v 		= 	S[1:(j-1),j]
		u 		= 	S[(j+1):m,j]
		mu[j]  	= 	mean(vcat(v[:],u[:]))
	end
	return 		mu
end

function trueordercanonicalform(	M;
									version		=	2,
									rev			=	false,
									firstval	=	1,
									factor		= 	false )
	if version 				== 	1
		# NB: this version only defined for direction == "up"
		m 					= 	length(M)
		n 					= 	length(unique(M))
		ocf 				= 	zeros(Int64,size(M)...)
		val 				= 	zeros(Float64,n)
		p 					= 	sortperm(M[:],alg=MergeSort)
		k 					= 	1
		val[k] 				= 	M[p[1]]
		ocf[p[1]] 			= 	k
		for i 				= 	2:m
			if 	M[p[i]]		> 	M[p[i-1]]
				k 				+=1
				val[k] 		= 	M[p[i]]
			end
			ocf[p[i]] 	= 	k
		end
	elseif version 			== 	2

		m 					= 	length(M)
		perm 				= 	sortperm(M[:],rev = rev,alg=MergeSort)
		oca 				= 	zeros(Int64,size(M)...)

		if m 				==	0
			return 				zeros(Int64,0),zeros(Int64,0)
		end

		if factor
			numvals 			= 	1
			post 				= 	1
			for p 				= 	1:m
				if M[perm[p]] 	!= 	M[perm[post]]
					post 		= 	p
					numvals 	+= 	1
				end
			end
			oca2rad 			= 	Array{Float64}(undef,numvals)
			oca2rad[1] 			= 	M[perm[1]]
		end

		post 					= 	[1]
		k 						= 	[1]
		trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
		if factor
			return oca,oca2rad
		else
			return oca
		end
	end
end

function trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
	for p 	 			= 	1:m
		if M[perm[p]] 	!= 	M[perm[post[1]]]
			post[1] 	= 	p
			k[1] 		=  	k[1]+1
			if factor
				oca2rad[k[1]] 	= 	M[perm[post[1]]]
			end
		end
		oca[perm[p]] 	= 	k[1]
	end
end



##########################################################################################

####	COMBINATIONS, PERMUTATIONS, AND SET OPERATIONS

##########################################################################################

function intervalcomplementuniqueunsortedinput(uniquepositiveintegers,intervalEndpoint)
	v = uniquepositiveintegers
	n = intervalEndpoint
	L = length(v)
	if L==0
		return 1:n
	elseif L==n
		return Array{Int64}(undef,0)
	else
		complementsupport = trues(n)
		complementsupport[v].=false
		complement = Array{Int64}(undef,n-L)
		marker = 0
		for i = 1:n
			if complementsupport[i]
				marker+=1
				complement[marker]=i
			end
		end
	end
	return complement
end

function integersinsameorder!(v::Array{Int64,1},maxradue::Int64)
	m = length(v)
	x = zeros(Int64,maxradue)
	for i = 1:m
		x[v[i]]+=1
	end
	y = Array{Int64}(undef,maxradue+1)
	y[1] = 1
	for i = 1:maxradue
		y[i+1]=y[i]+x[i]
	end
	for i = 1:length(v)
		u = v[i]
		v[i] = y[u]
		y[u]+=1
	end
	return v
end

function integersinsameorder(v::Array{Int64,1})
	# Returns the permutation z on {1,...,length(v)} such z[i]<z[j] iff either
	# (a) v[i] < v[j], or
	# (b) v[i] = v[j] and i < j
	if isempty(v)
		z = Array{Int64}(undef,0)
		return z
	else
		m = length(v)
		maxv = maximum(v)
		minv = minimum(v)
		minv = minv-1;
		x = zeros(Int64,maxv-minv)
		z = Array{Int64}(undef,length(v))
		for i = 1:m
			x[v[i]-minv]+=1
		end
		prevsum = 1
		for i = 1:length(x)
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = 1:m
			u = v[i]
			z[i] = x[u-minv]
			x[u-minv]+=1
		end
		return z
	end
end

function integersinsameorder!(v::Array{Int64,1})
	# Replaces v with the permutation z on {1,...,length(v)} such that z[i]<z[j] iff either
	# (a) v[i] < v[j], or
	# (b) v[i] = v[j] and i < j
	if isempty(v)
		z = Array{Int64}(undef,0)
		return z
	else
		m = length(v)
		maxv = maximum(v)
		minv = minimum(v)
		minv = minv-1;
		x = zeros(Int64,maxv-minv)
		z = Array{Int64}(undef,length(v))
		for i = 1:m
			x[v[i]-minv]+=1
		end
		prevsum = 1
		for i = 1:length(x)
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = 1:m
			u = v[i]
			v[i] = x[u-minv]
			x[u-minv]+=1
		end
	end
end



#=
- 	In beta; should be compared with integersinsameorderbycolumn3.  See
	/Users/greghenselman/Google Drive/GregDirectory/julia_gd/Julia/workshop/workshop_Oct2017.jl
-   Functionally equivalent to integersinsameorderbycolumn; returns a
	permutation z on {1,...,length(v)} so that for all j
	- cran(colptr,j) maps to cran(colptr,j), and
	- crows(colptr,v[z],j) is an array in sorted order
=#
function integersinsameorderbycolumn2(v::Array{Int64,1},colptr)
	numcols = length(colptr)-1
	m = length(v)
	v = v.-(minimum(v)-1)
	x = zeros(Int64,maximum(v))
	z = Array{Int64}(undef,length(v))
	for j = 1:numcols
		if colptr[j] == colptr[j+1]
			continue
		end
		for i = colptr[j]:(colptr[j+1]-1)
			x[v[i]]+=1
		end
		maxv = v[colptr[j]];   minv = maxv
		for i = (colptr[j]+1):(colptr[j+1]-1)
			if v[i] > maxv
			   maxv = v[i]
			elseif v[i] < minv
			   minv = v[i]
			end
		end
		prevsum = colptr[j]
		for i = minv:maxv
			sum = prevsum + x[i]
			x[i] = prevsum
			prevsum = sum
		end
		for i = colptr[j]:(colptr[j+1]-1)
			u = v[i]
			z[i] = x[u]
			x[u]+=1
		end
		for i = minv:maxv
			x[i] = 0
		end
	end
	return z
end



##########################################################################################

####	SEARCH SUBROUTINES

##########################################################################################

function getPairsLightWrite2!(
	rowval::Array{Tv,1},
	colptr::Array{Tv,1},
	rowfilt::Array{Tv,1},
	colfilt::Array{Tv,1},
	m::Integer,
	n::Integer,
	prows::Array{Tv,1},
	pcols::Array{Tv,1},
	numpairs::Array{Tv,1};
	verbose = false)  where Tv<:Integer

	col2firstplace = zeros(Tv,n)
	rowwisesum = zeros(Tv,m)

	if verbose
		println("starting search for pairs")
	end
	for j = 1:n
		firstplace = colptr[j]
		if firstplace < colptr[j+1]
			firstrow = rowval[firstplace]
			rowwisesum[firstrow]+=1
			if firstplace < colptr[j+1]-1
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newfilt = rowfilt[newrow]
					rowwisesum[newrow]+=1
					if newfilt > filt
						filt = newfilt
						firstplace = newplace
					end
				end
			end
			col2firstplace[j]=firstplace
		end
	end
	colwisesum = colsupportsum(colptr,n) # note allow extra on end
	colfiltptr = getcolptr2(colfilt,n)	# note allow extra on end
	colwisesumlinearized = integersinsameorderbycolumn2(colwisesum,colfiltptr)
	colnamesinorder = Array{Tv}(undef,n)
	colnamesinorder[colwisesumlinearized]=1:n
	ncoveredsupp = trues(m)
	pairmarker = 0
	for jp = 1:n
		j = colnamesinorder[jp]
		if col2firstplace[j]>0
			firstplace = col2firstplace[j]
			firstrow = rowval[firstplace]
			if firstplace==colptr[j+1]-1
				if ncoveredsupp[firstrow]	#######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			else
				firstweight = rowwisesum[firstrow]
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newweight = rowwisesum[newrow]
					if rowfilt[newrow]==filt && newweight <= firstweight && !ncoveredsupp[firstrow] && ncoveredsupp[newrow]
						firstweight = newweight
						firstrow = newrow
					end
				end
				if ncoveredsupp[firstrow] ######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			end
			for ip = cran(colptr,j)
				ncoveredsupp[rowval[ip]]=false
			end
		end
	end

	numpairs[1]=pairmarker
	if verbose
		print("done finding pairs")
	end
end

function finddownstreamelements_embeddedupperunitriangularmatrix(
	Mrv,
	Mcp,
	Mm,
	initialelements::Array{Int64,1},
	prows::Array{Int64,1},
	pcols::Array{Int64,1};
	verbose=false)

	if verbose
		print("starting to get downstream elements")
	end
	if length(prows) != length(pcols)
		print("length of p doesn't match length of q")
		return
	elseif length(prows)==0
		return Array{Int64}(undef,0)
	end
	n = length(prows)
	rowtranslator = Array{Int64}(undef,Mm)
	for i = 1:n
		rowtranslator[prows[i]]=i
	end
	prowsupp = falses(Mm)
	prowsupp[prows].=true
	downstreamsupport = falses(n)
	for i = 1:length(initialelements)
		row = initialelements[i]
		if prowsupp[row]
			downstreamsupport[rowtranslator[row]] = true
		end
	end
	for jp = n:-1:1
		j = pcols[jp]
		ran = cran(Mcp,j)
		for ip in ran
			rawrow = Mrv[ip]
			if prowsupp[rawrow] && downstreamsupport[rowtranslator[rawrow]]
				for kp in ran
					rawrow = Mrv[kp]
					if prowsupp[rawrow]
						downstreamsupport[rowtranslator[rawrow]] = true
					end
				end
				break
			end
		end
	end
	counter = 0
	for i = 1:n
		if downstreamsupport[i]
			counter+=1
		end
	end
	downstreamelements = Array{Int64}(undef,counter)
	counter = 0
	for i = 1:n
		if downstreamsupport[i]
			counter+=1
			downstreamelements[counter] = i
		end
	end
	if verbose
		print("done with downstream elements")
	end
	return downstreamelements
end

##########################################################################################

####	BARCODE UTILITIES

##########################################################################################


function nnzbars(D::Dict;dim = 1)
	sd = dim+2
	if !(0<sd<=length(D["plo"]))
		print("error: requested dimension is outside the range of calculated bars")
		return
	end
	plo = D["plo"][sd]
	phi = D["phi"][sd]
	lowfilt = D["grain"][sd-1]
	higfilt = D["grain"][sd]
	counter::Int64 = 0
	for i = 1:length(plo)
		if higfilt[phi[i]]!=lowfilt[plo[i]]
			counter+=1
		end
	end
	counter += length(D["tid"][sd])-length(plo)
	return counter
end


function barname2cyclename(D::Dict,barnumber = [1];dim = 1)
	if typeof(barnumber) <: Array
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		nummortals = length(plo)
		nzcycles = findall(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
		append!(nzcycles,Array((nummortals+1):length(tid)))
		return nzcycles[barnumber]
	elseif typeof(barnumber) <: UnitRange
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		nummortals = length(plo)
		nzcycles = findall(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
		append!(nzcycles,nummortals+1:length(tid))
		return nzcycles[barnumber]
	elseif typeof(barnumber)<:Number
		sd = dim+2
		tid = D["tid"][sd]
		plo = D["plo"][sd]
		phi = D["phi"][sd]
		numclasses = length(tid)
		nummortals = length(plo)
		counter = 0
		cyclename = 0
		for i = 1:nummortals
			if D["grain"][sd][phi[i]] != D["grain"][sd-1][plo[i]]
				counter+=1
			end
			if counter == barnumber
				cyclename = i
				break
			end
		end
		if cyclename == 0
			for i = (nummortals+1):length(tid)
				counter+=1
				if counter == barnumber
					cyclename = i
					break
				end
			end
		end
		return cyclename
	end
end

function getbetticurve(D::Dict,sd;ocf = false)
	# NB: sd = (homology dimension)+2
	# NB: D["ocg2rad"]: {grains > 0} --> rad
	# NB: the ocf barcode takes values in [0, #{grains}-1]

	if complexrank(D,dim=sd-2) ==0
		return Array{Float64}(undef,0,2)
	end

	# the plus 1 is for the zero grain, the "s" in ngrains is for "shifted",
	# since we shift up to avoid zero indices
	ngrains = length(D["ocg2rad"])
	v = zeros(Int64,ngrains)

	bco = barcode(D;dim = sd-2,ocf=true)
	# bco[bco.==Inf] = ngrains
	bco = convert(Array{Int64},bco)

	for i = 1:size(bco,1)
		#try
		ran = 1 .+ (bco[i,1]:(bco[i,2]-1))
		v[ran] = v[ran] .+ 1
		#catch e
		#	println(v[ran])
		#	println(bco[i,1]:(bco[i,2]-1))
		#	error(e)
		#end
	end

	if ocf == false
		u = sort(D["ocg2rad"])
	else
		u = Array(1:ngrains)
	end
	return hcat(u,v)
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
	append!(summands,[tid[sd][cyclenumber]])
	brv = ff2aflight(farfaces,firstv,sd-1,summands)
	supp = falses(length(farfaces[sd-2]))
	for k in brv
		supp[k] = !supp[k]
	end

	brv = findall(supp[tid[sd-1]])
	bcp = [1,length(brv)+1]
	brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
	brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

	plow2phigtranslator = zeros(Int64,numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	brv		= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
	brv 	= brv[findall(!iszero,brv)]

	brv = append!(brv,summands)

	return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
	if sd == 2
		# this case must be treated separately b/c julia indexing starts at 1
		numlowlows = 0
		numnlpll = 0
	else
		numlowlows = length(cp[sd-2])-1
		numnlpll = numlowlows-length(plo[sd-2])
	end
	numnlpl = length(cp[sd-1])-1-length(plo[sd-1])	# the -1 accounts for the fact that length(cp[sd-1]) = (# cells of dimension secard-1) - 1

	summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
	append!(summands,[tid[sd][cyclenumber]])

	if sd == 2
		return summands
	end

	supp = falses(numlowlows)
	sc = sd-1
	for j in summands	# this is a bit ridiculus; it's here bc i don't have a verified form of spmmF2 on hand
		for k in cran(cp[sc],j)
			i = rv[sc][k]
			supp[i] = !supp[i]
		end
	end
	brv 				= findall(supp[tid[sd-1]])
	bcp 				= [1,length(brv)+1]
	brv,bcp 			= spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpll)
	brv,bcp 			= spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpll)
	plow2phigtranslator = zeros(Int64,numlowlows)
	plow2phigtranslator[plo[sd-1]]=phi[sd-1]
	brv					= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
	brv 				= append!(brv,summands)
	return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Array{Int64,1})
	numclasses 	= length(cyclenumber)
	rep 	 	= Array{Array{Int64,1},1}(undef,numclasses)
	for p 		= 1:numclasses
		rep[p] 	= getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber[p])
	end
	return rep
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	if sd == 2
		numlowlows = 0
	else
		numlowlows = length(farfaces[sd-2])
	end
	numlows    = length(farfaces[sd-1])
	numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

	numclasses = length(cyclenumber)
	summands = Array{Array{Int64,1},1}(undef,numclasses)
	rep 	 = Array{Array{Int64,1},1}(undef,numclasses)
	summandsupp = falses(numlows)
	for i = 1:numclasses
		summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
		append!(summands[i],[tid[sd][cyclenumber[i]]])
		summandsupp[summands[i]].=true
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

		supp[:] .= false
		for j = 1:length(summands[i])
			for k = 1:m
				kk = lowfacemat[k,translator[summands[i][j]]]
				supp[kk] = !supp[kk]
			end
		end

		brv = findall(supp[tid[sd-1]])
		bcp = [1,length(brv)+1]
		brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numlowlows)
		brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numlowlows)

		rep[i] = append!(plow2phigtranslator[tid[sd-1][brv]],summands[i])
	end
	return rep
end

function getcycle(D::Dict,sd,cyclenumber)
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
	rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	return rrv
end

function getcycle(D::Dict,cyclenumber;dim = 1)
	if !haskey(D,"Lirv")
		println("Error: the function <getcycle(D::Dict,cyclenumber;dim = 1)> assumes a key value for \"Lirv\" in the input object D.  This key value is absent.")
		return
	end
	sd = dim+2
	Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
	Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];
	plo=D["plo"];phi=D["phi"];tid=D["tid"]
	if haskey(D,"farfaces")
		farfaces = D["farfaces"];firstv = D["firstv"];
		rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	else
		rrv = getcycle_cell(D["rv"],D["cp"],Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
	end
	return rrv
end

function numcols(cp)
	return length(cp)-1
end

function complexrank(C;dim=1)
	sd 		= 	dim+1
	if 	dim > C["input"]["maxdim"]+1 || dim < 0
		return 0
	elseif C["input"]["model"] == "complex"
		return length(C["cp"][sd])-1
	else
		return length(C["farfaces"][sd])
	end
end

function empteval(f,a,c)
	if isempty(a)
		return c
	else
		return f(a)
	end
end

##########################################################################################

####	USER-FRIENDLY UTILITIES

##########################################################################################


function boundarymatrices(C)
	if haskey(C,"farfaces")
		rv,cp = ff2complex(C["farfaces"],C["firstv"])
	else
		rv = C["rv"]
		cp = C["cp"]
	end
	return rv,cp
end

function classrep(
	D::Dict;
	dim = 1,
	class = 1,
	format = "vertex x simplex")

	if any(class.>nnzbars(D,dim=dim))
		print("error: the value for keyword argument <class> has an integer greater than the number of nonzero bars in the specified dimension")
		return
	elseif !(0<=dim<=D["input"]["maxdim"])
		print("error: barcodes were not computed in the specified dimension")
		return
	end
	if !haskey(D,"farfaces")
		format = "index"
	end

	if format == "vertex x simplex"
		return classrep_faces(D,dim = dim,class = class)
	elseif format == "vertex"
		return classrep_vertices(D,dim = dim, class = class)
	elseif format == "index"
		return classrep_cells(D,dim=dim,class=class)
	end
end

function classrep_cells(
	D::Dict;
	dim = 1,
	class = 1)

	sd = dim+2

	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		cyclename = barname2cyclename(D,class;dim=dim)
		rep = getcycle(D,sd,class)
	end

	return rep
end

function classrep_faces(
	D::Dict;
	dim = 1,
	class = 1)

	sd 				= 	dim+2
	rep 			= 	classrep_cells(D,dim=dim,class=class)

	vrealization 	= 	vertexrealization(D::Dict,sd-1,rep)
	vrealization 	= 	D["nvl2ovl"][vrealization]
	return vrealization
end

function classrep_vertices(
	D::Dict;
	dim = 1,
	class = 1)

	sd 			= 	dim+2
	rep 		= 	classrep_cells(D,dim=dim,class=class)

	vertices 	= 	incidentverts(D::Dict,sd-1,rep)
	vertices 	= 	D["nvl2ovl"][vertices]
	return vertices
end


function printval(var,varname)
	println(string(varname," = $(var)"))
end

function barcode(D::Dict;dim = 1,ocf = false)
	if 	haskey(D,:perseusjlversion)
		return	barcode_perseus(D,dim=dim)
	elseif haskey(D,:barcodes)
		return  D[:barcodes][dim+1]
	elseif !haskey(D,"cp") & !haskey(D,"farfaces")
		print("unrecognized object:")
		display(D)
		return
	elseif dim > D["input"]["maxdim"]
		maxdim 	= 	D["input"]["maxdim"]
		println("error: barcodes were not computed in dimensions greater than $(maxdim).")
		return
	end
	sd = dim+2
	plo = D["plo"][sd]
	phi = D["phi"][sd]
	tid = D["tid"][sd]
	lg 	= D["grain"][sd-1]
	hg  = D["grain"][sd]

	mortalprimagrain 	= 	lg[plo]
	mortalultragrain 	= 	hg[phi]

	finind 				= 	findall(mortalprimagrain .!= mortalultragrain)
	numfin 				= 	length(finind)
	numinf 				= 	length(tid)-length(plo)

	mortalprimagrain 	= 	mortalprimagrain[finind]
	mortalultragrain 	= 	mortalultragrain[finind]

	mortalran 			= 	1:numfin
	evergrran 			= 	numfin+1:numfin+numinf
	finran 				= 	1:(2*numfin+numinf)
							# stands for finite range; these are the linear
							# indices of the barcode array that take finite
							# values
	evrgrbran 			= 	length(tid)-numinf+1:length(tid)
							# stands for evergreen birth range; this satisfies
							# tid[ebergrbran] = [array of evergreen classes]

	bc 					= 	zeros(Int64,numfin+numinf,2)
	bc[mortalran,1]    .= 	mortalprimagrain
	bc[mortalran,2]    .= 	mortalultragrain
	bc[evergrran,1] 	= 	lg[tid[evrgrbran]]

	if !ocf
		bcc 				= 	copy(bc)
		bc 					= 	Array{Float64}(bc)
		bc[finran] 			= 	D["ocg2rad"][bcc[finran]]
		bc[evergrran,2]    .= 	Inf
	else
		bc 					= 	length(D["ocg2rad"]).-bc
	end

	return bc
end

function getpersistencediagramprimitives(
	C;
	dim = 1,
	ocf = false,
	descriptivetext = true,
	showsize = true)

	if haskey(C,"cyclerep")
		showsize = true
	else
		showsize = false
	end

	bco = barcode(C,dim = dim,ocf = ocf)
	rows = unique(bco,dims=1)
	numrows = size(rows,1)
	numbrs = size(bco,1)
	sd = dim+2

	if numbrs == 0
		x0=[];y0=[];l0=[];x1=[];y1=[];l1=[];x2=[];y2=[]
		return x0,y0,l0,x1,y1,l1,x2,y2
	end

	if showsize
		barsizes = Array{Int64}(undef,numbrs)
		for i = 1:numbrs
			barsizes[i] = length(C["cyclerep"][sd][i])
		end
	end

	if descriptivetext
		if showsize
			labels = fill("class/size  ",numrows)
		else
			labels = fill("class  ",numrows)
		end
	else
		labels = fill("",numrows)
	end

	D = Dict()
	for i = 1:numrows
		D[Symbol(rows[i,:])] = i
	end
	for j = 1:numbrs
		uniquerowind = D[Symbol(bco[j,:])]
		if labels[uniquerowind] in ["", "class  ","class/size  "]
			if showsize
				labels[uniquerowind] = string(labels[uniquerowind],"$(j)/$(barsizes[j])")
			else
				labels[uniquerowind] = string(labels[uniquerowind],"$(j)")
			end
		else
			if showsize
				labels[uniquerowind] = string(labels[uniquerowind],", $(j)/$(barsizes[j])")
			else
				labels[uniquerowind] = string(labels[uniquerowind],", $(j)")
			end
		end
	end
	if any(rows[:,2].!=Inf)
		topheight = 1.1*maximum(rows[:,2][rows[:,2].!=Inf])
	else
		topheight = 0
	end
	infrows = findall(rows[:,2].==Inf)
	finrows = findall(rows[:,2].!=Inf)

	x0 = rows[infrows,1]
	y0 = x0
	l0 = labels[infrows]
	x1 = rows[finrows,1]
	y1 = rows[finrows,2]
	l1 = labels[finrows]
	x2 = [minimum(rows[:,1]),maximum(rows[:,1])]
	y2 = [topheight,topheight]
	return x0,y0,l0,x1,y1,l1,x2,y2
end

function plotpersistencediagram_pjs(C;dim = 1,showlabels = false,ocf = false)
	if showlabels == true
		mmode = "markers+text"
	else
		mmode = "markers"
	end

	x0,y0,l0,x1,y1,l1,x2,y2 = getpersistencediagramprimitives(
		C;
		dim = dim,
		ocf = ocf)

	if length(x0)==0 && length(x1) == 0
		print("There are no persistent homology classes in dimension $(dim).")
		return
	end

    trace0 = PlotlyJS.scatter(
    	x=x0,
    	y=y0,
        text=l0,
        mode=mmode,
        marker_color = "red",
        textposition="top center",
        hoverinfo = "x+text",
        marker_size=6,
        textfont_family="Raleway, sans-serif")
    trace1 = PlotlyJS.scatter(
    	x=x1,
    	y=y1,
        text=l1,
        mode=mmode,
        marker_color = "black",
        textposition="top center",
        hoverinfo = "x+y+text",
        marker_size=5,
        textfont_family="Raleway, sans-serif")

    L = makelayout_pjs(trace1,showlegend = false)
    PlotlyJS.plot([trace1,trace0],L)
end

function classrep_pjs(
	D::Dict;
	dim = 1,
	class=1,
	showcloud= [],
	coords = [],
	model= D["input"]["model"],
	embeddingdim = 3,
	embeddingobj = [],
	specrange = [-Inf,Inf],
	classcolor = "spectral",
	cloudcolor = [],
	textlabels = [],
	showlabels = false,
	showedges = true,
	alwaysshowcyclelabels = false)

	###
	if !(embeddingdim in [2,3])
		println("Invalid key value: <embeddingdim> must be either 2 or 3.")
		return
	end
	if embeddingobj == "hop" && showcloud == true
		println("Error: keyword <embeddingobj> may only take value <\"hop\"> when keyword <showcloud> takes value <false>.")
		println()
		return
	end

	###
	sd 	= dim+2
	facecard 	= dim+1

	##
	if showcloud == []
		showcloud = true
	end
	if embeddingobj == []
		embeddingobj = "dmat"
	end

	###
	cyclename = barname2cyclename(D,class;dim = dim)
	if haskey(D,"cyclerep")
		rep = D["cyclerep"][sd][class]
	else
		rep = getcycle(D,sd,cyclename)
	end
	classvinnewspace = incidentverts(D,sd-1,rep)
	classvinoldspace = D["nvl2ovl"][classvinnewspace]
	if showcloud == true
        # subselect farfaces?
		vsupp = trues(length(D["farfaces"][1]))
		vsupp[classvinnewspace] .= false
		compvinnewspace = findall(vsupp)
		compvinoldspace = D["nvl2ovl"][compvinnewspace]
	else
		compvinoldspace = []
	end

	###
	if haskey(D["input"],"pointlabels")
		textlabels = D["input"]["pointlabels"]
	else
		m = length(D["farfaces"][1])
		textlabels = Array{String,1}(m)
		for i = 1:m
			textlabels[i] = "$(i)"
		end
	end
	if showlabels == true || showlabels == "cycle" || showlabels == "all"
		showlabtemp = true
	else
		showlabtemp = false
	end

	###
	if (classcolor == "spectral") | (embeddingobj == "hop")
		vrealization = vertexrealization(D,dim=dim,class=class)
 		vrealization = D["nvl2ovl"][vrealization]
		vertexinverter = Array{Int64}(undef,maximum(classvinoldspace))
		vertexinverter[classvinoldspace]=1:length(classvinoldspace)
		classedges = d1faces(vrealization)
		edges_orderverts = vertexinverter[classedges]
		L = graphlaplacian_normal(edges_orderverts)
		efact_class = eigen(L)
		efact_class = efact_class.vectors[:,1:4]
	end

	if showcloud && embeddingobj == "hop"
		hopedges = findall(hoprange[1]<=D["grain"][2] & D["grain"][2].<=hoprange[2])
		cloudedges = vertexrealization(D,dim=1,hopedges)
		cloudedges_orderverts = vetexinverter[cloudedges]
	end

	if coords == []
		if D["input"]["pc"] == "n/a"
			print("No point cloud is available.  Please consider using the mds keyword argument to generate a Euclidean embedding from the distance matrix (see documentation).")
			return "nopointcloud","nopointcloud"
		elseif D["input"]["model"]=="pc"
			coords = D["input"]["genera"]
		elseif typeof(D["input"]["pc"]) <: Array{Float64} || typeof(D["input"]["pc"]) <: Array{Int64}
			coords = D["input"]["genera"]
		end
		if size(coords,1) > 3
			print("The input coordinates have dimension greater than 3.  The generated plot will use the first three to represent each point. For more options re: high-dimensional representations, please see the documentation for multidimensional scaling.")
			coords = coords[1:3,:]
		end
		if !showcloud
			coords = coords[:,classvinoldspace]
		end
	elseif coords == "mds"
		if showcloud
			if embeddingobj == "dmat"
				metricmatrix = D["input"]["genera"]
				metricmatrix = metricmatrix .- minimum(metricmatrix)
				for i = 1:size(metricmatrix,1)
					metricmatrix[i,i]=0
				end
			elseif embeddingobj == "hop"
				metricmatrix = hopdistance(cloudedges_orderverts,inputis = "edges")
			end
		else
			if embeddingobj == "dmat"
				metricmatrix = D["input"]["genera"][classvinoldspace,classvinoldspace]
				metricmatrix = metricmatrix - minimum(metricmatrix)
				for i = 1:size(metricmatrix,1)
					metricmatrix[i,i]=0
				end
			elseif embeddingobj == "hop"
				metricmatrix = hopdistance(edges_orderverts,inputis = "edges")
			end
		end
		coords = transform(fit(MDS, float.(metricmatrix), maxoutdim = embeddingdim, distances=true)) #classical_mds(metricmatrix,embeddingdim)
		# coords = round.(coords,10)
		model = "pc"
	end

	###
	if !showcloud
		textlabels = textlabels[classvinoldspace]
		classvinoldspace = 1:length(classvinoldspace)
	end

	###
	if classcolor == "spectral"
		classcolor1 = efact_class[:,2]
		classcolor1 = (classcolor1 .-minimum(classcolor1))  ./  (maximum(classcolor1)-minimum(classcolor1))
		classcolor2 = efact_class[:,3]
		classcolor2 = (classcolor2 .-minimum(classcolor2))  ./  (maximum(classcolor2)-minimum(classcolor2))
	end

	###
	T1 = maketrace_pjs(
		coords;
		model = model,
		subset = classvinoldspace,
		textlabels = textlabels,
		showlabels = showlabtemp)
	T1["marker_cmin"] = 0
	T1["marker_cmax"] = 1
	T1["marker_color"] = classcolor1
	T1["marker_line_color"] = classcolor2
	T1["marker_line_width"] = 2
	T1["marker_colorscale"] = "Jet"
	T1["marker_line_colorscale"] = "Jet"
	data = [T1]

	if !isempty(compvinoldspace)
		if showlabels == "all"
			showlabtemp = true
		else
			showlabtemp = false
		end
		T2 = maketrace_pjs(
			coords;
			model = model,
			subset = compvinoldspace,
			color = "rgb(31,119,180)",
			opacity = 0.5,
			textlabels = textlabels,
			showlabels = showlabtemp)
		append!(data,[T2])
	end

	if showedges
		faces = classrep(D,dim=dim,class=class)
		edges = d1faces(faces)
		if showcloud
			T3 =  edgetrace_pjs(coords,edges,model=model	)
		else
			# by default vertices in the edge list index into the original
			# point cloud; if we're not showing the ambient cloud, then we have
			# by this point deleted the "irrelevant" columns from the coordinate
			# matrix, which changes the indices.  padding with zeros corrects
			# for this difference
			# note that for design reasons we re-defined classvinoldspace above,
			# so must re-define it
			classvinoldspace_original = D["nvl2ovl"][classvinnewspace]
			coords_padded = zeros(size(coords,1),maximum(classvinoldspace_original))
			coords_padded[:,classvinoldspace_original] = coords
			T3 =  edgetrace_pjs(coords_padded,edges,model=model	)
		end

		append!(data,T3)
	end

	if model == "pc"
		dim = size(coords,1)
	else
		dim = size(coords,2)
	end
	xspan = maximum(T1["x"])-minimum(T1["x"])
	yspan = maximum(T1["y"])-minimum(T1["y"])
	xa = xspan/xspan
	ya = yspan/xspan
	if dim == 3
		zspan = maximum(T1["z"])-minimum(T1["z"])
		za = zspan/xspan
	end

	layout = makelayout_pjs(data)

	return data, layout

end

function plotclassrep_pjs(
	D::Dict;
	dim = 1,
	class=1,
	showcloud = [],
	coords = [],
	model= "pc",
	embeddingdim = 3,
	embeddingobj = [],
	specrange = [-Inf,Inf],
	classcolor = "spectral",
	cloudcolor = [],
	textlabels = [],
	showlabels = "cycle",
	showedges = false)

	if D["input"]["model"]!="pc" && D["input"]["pc"] == "n/a" && coords == []
		print("No point cloud is available.  Coordinates may be supplied by the user with the <coords> keyword argument, or generated automatically via the mds keyword argument.  Please see documentation.")
		return	"nopointcloud","nopointcloud"
	end

	data,layout = classrep_pjs(
		D;
		dim = dim,
		class=class,
		showcloud = showcloud,
		coords = coords,
		model="pc",
		embeddingdim = embeddingdim,
		embeddingobj = embeddingobj,
		specrange = specrange,
		classcolor = classcolor,
		cloudcolor = cloudcolor,
		textlabels = textlabels,
		showlabels = showlabels)

	if data != "nopointcloud"
		return PlotlyJS.plot(data,layout)
	else
		return
	end
end

function betticurve(D::Dict;dim = 1,ocf = false)
	sd = dim+2
	return getbetticurve(D,sd,ocf = ocf)
end

function plotbetticurve_pjs(D::Dict;dim=1,ocf = false)
	bcu = betticurve(D;dim = dim, ocf = ocf)
	T = PlotlyJS.scatter(x = bcu[:,1],y=bcu[:,2],mode = "line",line_width=1)
	L = makelayout_pjs(T)
	PlotlyJS.plot(T,L)
end

##########################################################################################

####	BARCODE PLOTTING
####	Contributed by Paul Breiding and Sara Kališnik Verovšek, 25 January 2018

##########################################################################################

"""
	plotbarcode_pjs(C; <keyword arguments>)

Generate barcode diagram in PlotlyJS for a dictionary object `C` returned by
the function `eirene`.

# Keyword Arguments
- `dim = 0:C["input"]["maxdim"]`: (homological) homological dimensions to plot
- `sortby = "birth"`: determines whether bars appear in order of `"birth"` (that is, birth time) or `"age"` (that is, death time - birth time)
- `minage = zeros(Float64,length(dim))`: array or range of real numbers (tolerances), one for each dimension to be plotted; bars that do not meet this minimum lifetime requirement will not be plotted
- `lw = 2`: line width of bars to be plotted
"""
function plotbarcode_pjs(C::Dict; dim = 0:C["input"]["maxdim"], sortby = "birth", minage = zeros(Float64,length(dim)), lw=2)

        @assert length(dim) == length(minage) "Number of dimensions (was $(length(dim))) must be the same as the number of minageerance values (was $(length(minage)))."

        # colors = Colors.distinguishable_colors(length(dim)+1,[RGB(1,1,1)])[2:end]
        cols = Colors.colormap("Blues", mid = 0.5)
        rangevar = Int.(round.(collect(range(50,stop=100,length=length(dim)+1))))
        colors = map(r -> cols[r], rangevar)

        B = [barcode(C, dim = d) for d in dim]
        #upper_limit = maximum([maximum(b[b.< Inf]) for b in B])
		upper_limit = [empteval(maximum,b[b.< Inf],0) for b in B]
		upper_limit = empteval(maximum,upper_limit,0)

        traces = map(1:length(B)) do j
             b = B[j]
             lengths = b[:,2]-b[:,1]
             b = b[lengths .> minage[j],:]
             if size(b,1) == 0
                 return [PlotlyJS.scatter(;x=[0,0], y=[0,0], mode="lines",  line_width = lw, line_color = colors[j], name = "dimension $(dim[j])")]
             end
             s = sortperm(b[:,2]-b[:,1],alg=MergeSort)
             b = b[s,:]

             i = findall(x->x==Inf, b[:,2])
             b[i,2] .= 2 * upper_limit

             if sortby == "age"
             elseif sortby == "birth"
                s = sortperm(b[:,1],alg=MergeSort)
                b = b[s,:]
             else
                println("The second argument must be either \"length\" or \"lowerlimit\".")
                return 0
             end

             return [PlotlyJS.scatter(;x=b[1,:], y=[1,1], mode="lines",  line_width = lw, line_color = colors[j], name = "Dimension $(dim[j])");
             [PlotlyJS.scatter(;x=b[i,:], y=[i,i], mode="lines",  line_width = lw, line_color = colors[j], showlegend=false) for i in 2:size(b,1)]]
         end

         for i = 2:length(traces)
             for j = 1:length(traces[i])
                 traces[i][j][:y] = traces[i][j][:y] + traces[i-1][end][:y] .+ 10
             end
         end
         traces = vcat(traces...)

         x = maximum(vcat([t[:x] for t in traces]...))
         y = maximum(vcat([t[:y] for t in traces]...))

         layout = PlotlyJS.Layout(;
         xaxis = attr(rangevar = [-.001, x+0.001], showgrid=false, zeroline =false, title = "ϵ"),
         yaxis = attr(rangevar = [0,y+0.1], showgrid=false, ticks = false))
         return PlotlyJS.plot(traces[end:-1:1], layout)
end

##########################################################################################

####	SPECTRAL FUNCTIONS

##########################################################################################

function submatrixsublevellaplacianeigenstats(A;indices=1:size(A,1),threshold = Inf)
	m = length(indices)
	L = A[indices,indices]
	for alpha in L
		if alpha < 0
			print("Please note, it appears the Laplacian is being generated from a matrix with one or more negative entries.")
		end
	end
	for j = 1:m
		for i = 1:m
			if L[i,j] > threshold
				L[i,j] = 0
			end
		end
	end
	for i = 1:m
		L[i,i]=0
	end
	c = sum(L,1)
	for i = 1:m
		if c[i] != 0
			c[i] = c[i]^(-1/2)
		end
	end
	for i = 1:m
		L[i,i] = -c[i]
	end
	L = -L
	L = broadcast(*,L,c)
	L = broadcast(*,L,c')
	L = Symmetric(L)

	F = eigen(L)
	return F
end

function graphlaplacian_normal(edges)
	numverts = maximum(edges)
	numedges = size(edges,2)

	d = zeros(Int64,numverts)
	for k in edges
		d[k]+=1
	end
	d = d.^(-1/2)
	L = zeros(numverts,numverts)
	for k = 1:numedges
		i = edges[1,k]
		j = edges[2,k]
		coeff = -d[i]*d[j]
		L[i,j]=coeff
		L[j,i]=coeff
	end
	for k = 1:numverts
		L[k,k]=1
	end
	return Symmetric(L)
end

function pcloudevec(pcloud;indices=1:size(pcloud,2),threshold = Inf,eval = 2)
	pcloud = pcloud[:,indices]

	l = length(indices)
	if isodd(l)
		A = zeros(l+1,l+1)
		A[1:l,1:l] = Distances.pairwise(Euclidean(),pcloud,dims=2)
		eval = 3
	else
		A = Distances.pairwise(Euclidean(),pcloud,dims=2)
	end

	F = submatrixsublevellaplacianeigenstats(A,threshold = threshold,statrange = eval:eval)
	v = F.vectors
	if isodd(l)
		v = v[1:end-1]
	end
	v = v[:]
	v = v-minimum(v)
	v = v/maximum(v)
	return v
end

##########################################################################################

####	GRAPHING

##########################################################################################

function maketrace_pjs(
	pcloud;
	model = "pc",
	subset = 1:size(pcloud,2),
	color = "rgb[250,250,250]",
	colorscale = "Jet",
	threshold = Inf,
	opacity = 1,
	outlineonly = false,
	textlabels = [],
	markersize = 5,
	showlabels = false)

	if model == "pc"
		dim = size(pcloud,1)
		x = pcloud[1,subset]
		y = pcloud[2,subset]
		if dim >= 3
			z = pcloud[3,subset]
			if dim > 3
				print("It appears the dimension of the input point cloud exceeds 3. Using the first three coordinates only.\n")
			end
		end
	elseif model == "points"
		dim = size(pcloud,2)
		x = pcloud[subset,1]
		y = pcloud[subset,2]
		if dim >= 3
			z = pcloud[subset,3]
			if dim > 3
				print("It appears the dimension of the input point cloud exceeds 3. Using the first three coordinates only.\n")
			end
		end
	end
	x = x[:]
	y = y[:]
	if dim == 3
		z = z[:]
	end

	if outlineonly
		symb = "circle-open"
	else
		symb = "circle"
	end

	if dim == 2
		T = PlotlyJS.scatter(x = x,y = y)
	else
		T = PlotlyJS.scatter3d(;x = x, y= y, z=z,marker_line_width=2)
	end

	T["marker_size"] = markersize
	T["marker_opacity"] = opacity
	T["marker_symbol"] = symb
	T["autocolorscale"] = false

	if !isempty(textlabels)
		T["text"] = textlabels[subset]
	elseif showlabels
		T["text"] = Array{String,1}(length(x))
		for i = 1:length(x)
			T["text"][i] = "$(i)"
		end
	end
	if showlabels
		T["mode"] = "markers+text"
	else
		T["mode"] = "markers"
	end

	if isa(color,String) && color != "spectral"
		T["marker_color"] = color
		T["marker_line_color"] = color
	else
		if color == "spectral"
			if model == "points"
				pcloud = pcloud'
			end
			v = pcloudevec(pcloud;indices=subset,threshold = threshold,eval = 2)
		elseif isa(color,Array)
			v = color
		end
		v = v - minimum(v)
		v = v/maximum(v)
		T["marker_cmin"] = 0
		T["marker_cmax"] = 1
		T["marker_color"] = v
		T["marker_colorscale"] = colorscale
	 	T["marker_line_color"] = "green"
	end
	return T
end

function makelayout_pjs(data;showlegend = false,equalaxis = true,scenecolor = "rgb(2,44,82)")
	if !(typeof(data)<:Array)
		T = data
		data = Array{Any}(undef,1)
		data[1] = T
	end

	if data[1]["model"] == "scatter3d"
		dim=3
	else
		dim=2
	end

	if equalaxis
		xmax = -Inf
		xmin = Inf
		ymax = -Inf
		ymin = Inf
		if dim == 3
			zmax = -Inf
			zmin = Inf
		end

		for i = 1:length(data)
			T = data[i]
			if !isempty(T["x"])
				xmax = max(xmax,maximum(T["x"]))
				xmin = min(xmin,minimum(T["x"]))
				ymax = max(ymax,maximum(T["y"]))
				ymin = min(ymin,minimum(T["y"]))
				if dim == 3
					zmax = max(zmax,maximum(T["z"]))
					zmin = min(zmin,minimum(T["z"]))
				end
			end
		end
		if xmin == Inf
			L = PlotlyJS.Layout()
			return L
		end
		xspan = xmax - xmin
		yspan = ymax - ymin
		xa = xspan/xspan
		ya = yspan/xspan
		if dim == 3
			zspan = zmax-zmin
			za = zspan/xspan
		end
	else
		xa=1
		ya=1
		if dim==3
			za=1
		end
	end

	if dim==3
		L = PlotlyJS.Layout(
			showlegend = showlegend,
			width=1000,
			height=800,
			margin=attr(l=50, r=50, b=50, t=50),
			scene = attr(
				aspectratio=attr(x=xa,y=ya,z=za),
				aspectmode = "manual"
			),
		)
	else
		L = PlotlyJS.Layout(
			showlegend = showlegend,
			hovermode = "closest",
			width=1000,
			height=800,
			margin=attr(l=50, r=50, b=50, t=50),
			scene = attr(
				aspectratio=attr(
					x=xa,y=ya
				),
				aspectmode = "manual"
			),
		)
	end

	return L
end

function ezplot_pjs(
	pcloud;
	model = "pc",
	subset = 1:size(pcloud,2),
	color = "rgb[250,250,250]",
	colorscale = "Jet",
	threshold = Inf,
	opacity = 1,
	outlineonly = true,
	textlabels = [],
	markersize = 5,
	showlabels = false)

	T = maketrace_pjs(
		pcloud;
		model = model,
		subset = subset,
		color = color,
		colorscale = colorscale,
		threshold = threshold,
		opacity = opacity,
		outlineonly = outlineonly,
		textlabels = textlabels,
		markersize = markersize,
		showlabels = showlabels)

	L = makelayout_pjs(T)

	PlotlyJS.plot(T,L)
end

function edgetrace_pjs(coordinates,edges;model="pc")
	if model != "pc"
		coordinates = coordinates'
	end

	edgetraces = []
	if size(coordinates,1) == 2
		for i = 1:size(edges,2)
			verts = edges[:,i]
			if edgetraces == []
				trace = PlotlyJS.scatter(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					line = attr(color="#1f77b4", width=1.5),
					mode = "lines",
					name = "Dim 1 Faces")
				edgetraces = [trace]
			else
				trace = PlotlyJS.scatter(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					line = attr(color="#1f77b4", width=1.5),
					mode = "lines",
					name = "edge ($(verts[1]),$(verts[2]))",
					showlegend = false)
				append!(edgetraces,[trace])
			end
		end
	elseif size(coordinates,1) == 3
		for i = 1:size(edges,2)
			verts = edges[:,i]
			if edgetraces == []
				trace = PlotlyJS.scatter3d(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					z = coordinates[3,verts],
					line=attr(color="#1f77b4", width=1.5),
					mode = "lines",
					opacity = 0.5,
					name = "Dim 1 Faces")
				edgetraces = [trace]
			else
				trace = PlotlyJS.scatter3d(
					x = coordinates[1,verts],
					y = coordinates[2,verts],
					z = coordinates[3,verts],
					line=attr(color="#1f77b4", width=1.5),
					name = "edge ($(verts[1]),$(verts[2]))",
					showlegend = false,
					mode = "lines",
					opacity = 0.5)
				append!(edgetraces,[trace])
			end
		end
	else
		print("It appears the coordinates provided have dimension other than 2 or 3; these are the only two currently supported.")
	end
	return edgetraces
end

function d1faces(facesbycol)
	facecard = size(facesbycol,1)
	numfaces = size(facesbycol,2)
	M = maximum(facesbycol)
	supp = falses(M)
	for m in facesbycol
		supp[m] = true
	end
	vertices = findall(supp)
	numverts = length(vertices)
	translator = Array{Int64}(undef,M)
	translator[vertices]=1:numverts
	supp = falses(numverts,numverts)
	for i = 1:facecard-1
		for j = (i+1):facecard
			for k = 1:numfaces
				row = translator[facesbycol[i,k]]
				col = translator[facesbycol[j,k]]
				supp[max(row,col),min(row,col)]=true
			end
		end
	end
	numedges = count(!iszero,supp)
	edges = Array{Int64}(undef,2,numedges)
	counter = 1
	for j = 1:numverts
		for i = (j+1):numverts
			if supp[i,j]
				edges[1,counter]=j
				edges[2,counter]=i
				counter+=1
			end
		end
	end
	return vertices[edges]
end

##########################################################################################

####	MISC

##########################################################################################

function pairwiseisequal(X;under=identity,returnarray=false)
	l 					= 	length(X)
	for p 				= 	1:l
		if !isequal(under(copy(X[1])),under(copy(X[p])))
			if returnarray
				A 		=	trues(length(X),length(X))
				for p   =	1:length(X)
					for q 	=	(p+1):length(X)
						A[p,q] 	=	isequal(under(copy(X[p])),under(copy(X[q])))
					end
				end
				return 		A
			else
				return 		false
			end
		end
	end
	return 					true
end

# stands for extension-by-constant
function ec(v,p,k)
	if 0 < p <= length(v)
		return v[p]
	else
		return k
	end
end

function ezread(s)
	if s[end-2:end] == "csv"
		return readdlm(s,',','\r')
	elseif s[end-2:end] == "txt"
		return readdlm(s)
	else
		println("Please ensure the input file is either comma separated (.csv)
		or space delimited (.prn)")
	end
end

function ezlabel(y)
	l = length(y)
	x = Array{String,1}(undef,l)
	for i = 1:l
		x[i] = try
			convert(String,y[i])
		catch
			"$(i)"
		end
	end
	return x
end

function separatelabels(s,side)
	if side == "left" || side == "right"
		numlabels = size(s,1)
	else
		numlabels = size(s,2)
	end
	labels = Array{String,1}(numlabels)

	if side == "none"
		for i = 1:numlabels
			labels[i] = "$i"
		end
	end

	if side == "left"
		for i = 1:numlabels
			labels[i] = "$(s[i,1])"
		end
		s = s[:,2:end]
	elseif side == "right"
		for i = 1:numlabels
			labels[i] = "$(s[i,end])"
		end
		s = s[:,1:end-1]
	elseif side == "top"
		for i = 1:numlabels
			labels[i] = "$(s[1,i])"
		end
		s = s[2:end,:]
	elseif side == "bottom"
		for i = 1:numlabels
			labels[i] = "$(s[end,i])"
		end
		s = s[1:end-1,:]
	end
	return s,labels
end

function csvreadrow(fp;row=1,rowtype=Float64)
	M 					=	CSV.read(fp,header=0,skipto=row,limit=1)
	if rowtype 			==	nothing
		return M
	else
		M 				=	convert(Matrix{rowtype},M)
		M 				=	vec(M)
		return M
	end
end

function binom(x,y)
	k = 1;
	for i = x:-1:(x-y+1)
		k = i*k
	end
	for i = 2:y
		k = k/i
	end
	k = convert(Int64,k)
	return k
end

function yafterx!(y::Array{Tv,1},x::Array{Tv,1}) where Tv<:Integer
	for i = 1:length(x)
		x[i] = y[x[i]]
	end
end

function yafterx(y::AbstractVector{Tv},x) where Tv
	z = Array{Tv}(undef,length(x))
	for i = 1:length(x)
		z[i] = y[x[i]]
	end
	return z
end

function full2ss(A)
# 	added 12/27/2017; ss stands for 'sparse support'
# 	input: a full array A
# 	output: the support of A, encoded in sparse column format
	m,n = size(A)
	rv  = findall(!iszero,A)
	rv  = mod.(rv-1,m)+1
	cp  = zeros(Int64,n+1)
	cp[1] = 1
	for p = 1:n
		cp[p+1] = cp[p]+count(!iszero,A[:,p])
	end
	return rv,cp
end


function colsupportsum(colptr,n::Integer)
	x = Array{Int64}(undef,n)
	@inbounds begin
	for i = 1:n
		x[i] = colptr[i+1]-colptr[i]
	end
	end
	return x
end

function rowsupportsum(Arv,Acp,Am::Int64,cols)
	x = zeros(Int64,Am)
	for j in cols
		for ip in cran(Acp,j)
			x[Arv[ip]]+=1
		end
	end
	return x
end

function getcolptr2(orderedpositiveintegerlist::Array{Tv,1},howfartolookbeforestopping::Tv) where Tv<:Integer
	#### please note: order can be ascending or descending
	v = orderedpositiveintegerlist
	if isempty(v)
		return []
	end
	colptr = Array{Int64}(undef,length(v)+1)
	colptr[1] = 1
	transitioncounter = 1
	currentvalue = v[1]
	for i = 2:howfartolookbeforestopping
		if v[i]!=currentvalue
			transitioncounter+=1
			colptr[transitioncounter]=i
			currentvalue = v[i]
		end
	end
	colptr[transitioncounter+1]=howfartolookbeforestopping+1
	deleteat!(colptr,(transitioncounter+2):(length(v)+1))
	return colptr
end

function addinteger!(v::Array{Tv,1},k::Int64) where Tv
	for i = 1:length(v)
		v[i]+=k
	end
end

function sparseadjacencymatrix(A;inputis = "adjacencymatrix")
	if inputis == "adjacencymatrix"
		if size(A,1) != size(A,2)
			print("Error: unless the <inputis> keywork argument has value 'edges', the input array must be square.")
			return
		end
		m = size(A,1)
		rv = Array{Int64}(undef,0)
		cp = zeros(Int64,m+1)
		cp[1] = 1
		for i = 1:m
			adjverts = findall(!iszero,A[:,i])
			append!(rv,adjverts)
			cp[i+1] = cp[i]+length(adjverts)
		end
		return rv, cp
	elseif inputis == "edges"
		if size(A,1) != 2
			print("Error: when the <inputis> keywork argument has value 'edges', the input array must have exactly two rows.")
		end
		m = maximum(A)
		adjmat = falses(m,m)
		for i = 1:size(A,2)
			adjmat[A[1,i],A[2,i]]=true
		end
		adjmat[findall(transpose(adjmat))] .= true
		return sparseadjacencymatrix(adjmat)
	end
end

function hopdistance_sparse(rv,cp)
	m = length(cp)-1
	H = zeros(Int64,m,m)
	for i = 1:m
		c = 0
		metnodes = falses(m)
		metnodes[i] = true
		fringenodes = falses(m)
		fringenodes[i] = true
		fringelist = [i]

		while !isempty(fringelist)
			c+=1
			for j in fringelist
				for k in crows(cp,rv,j)
					if !metnodes[k]
						metnodes[k] = true
						fringenodes[k] = true
						H[k,i] = c
					end
				end
			end
			fringelist = findall(fringenodes)
			fringenodes[:].= false
		end
		H[.!metnodes,i].=m+1
	end
	return H
end

function hopdistance(rv,cp)
	return hopdistance_sparse(rv,cp)
end

function hopdistance(A;inputis = "fulladj")
	if inputis == "fulladj"
		rv,cp = sparseadjacencymatrix(A)
	elseif inputis == "edges"
		rv,cp = sparseadjacencymatrix(A,inputis="edges")
	else
		rv,cp = A
	end
	return hopdistance_sparse(rv,cp)
end


##########################################################################################

####	VIETORIS-RIPS CONSTRUCTION

##########################################################################################

function buildcomplex3(symmat::Array{Tv},maxsd; dictionaryoutput = true, verbose = false) where Tv

	grain = Array{Array{Int64,1}}(undef,maxsd+1)
	farfaces = Array{Array{Int64,1}}(undef,maxsd+1)
	prepairs = Array{Array{Int64,1}}(undef,maxsd+1)
	firstv = Array{Array{Int64,1}}(undef,maxsd+1)

	farfaces[maxsd+1] = Array{Int64}(undef,0)
	firstv[maxsd+1] = ones(Int64,1)
	grain[maxsd+1] = Array{Int64}(undef,0)
	prepairs[maxsd+1] = Array{Int64}(undef,0)

	m = size(symmat,1)
	w = vec(offdiagmean(symmat,defaultvalue=0)) 	# modified 02/12/2018

	vperm = sortperm(-w,alg=MergeSort)
	symmat = symmat[vperm,vperm]

	farfaces[1] = convert(Array,1:m)
	firstv[1] = convert(Array,1:(m+1))
	grain[1] = diag(symmat)
	prepairs[1] = Array{Int64}(undef,0)

	r,c,z = generate2faces(symmat)
	farfaces[2] = r
	firstv[2] = c
	grain[2] = z
	prepairs[2] = Array{Int64}(undef,0)

	if maxsd == 3
		generate3faces!(farfaces,firstv,grain,prepairs,m,symmat;verbose = verbose)
		if dictionaryoutput == true
			D = Dict{String,Any}(
				"farfaces" => farfaces,
				"firstv" => firstv,
				"grain" => grain,
				"prepairs" => prepairs,
				"symmat" => symmat,
				"nvl2ovl"=>vperm)
			return D
		else
			return farfaces,firstv,grain,prepairs,symmat,vperm
		end
	end

	fpi = Array{Int64}(undef,0)
	ff2pv = Array{Int64}(undef,0)
	pmhist = zeros(Int64,m,m)

	for sd = 3:maxsd
		if verbose
			print(["set cardinality = " sd])
			println(["num sd-1 cells" length(farfaces[sd-1])])
		end

		nl = length(farfaces[sd-1])
		nll = length(farfaces[sd-2])

		startlength = nl
		stepsize = min(10^7,Int(ceil(nl/4)))

		npsupp = trues(nl)
		pflist = Array{Int64}(undef,nl)
		jrv = farfaces[sd-1]
		jcp = firstv[sd-1]
		jz = grain[sd-1]
		zll= grain[sd-2]
		izfull = Array{Int}(undef,nll)
		r = Array{Int64}(undef,startlength)
		z = Array{Int64}(undef,startlength)
		c = Array{Int64}(undef,m+1)
		c[1]=1
		numpairs = [0]
		facecount = [0]
		if sd == maxsd-1
			ff2pv = Array{Int64}(undef,nl)
			ff2pv[:].= m+1
		end
		if sd == maxsd
			#### sort j-matrix by grain
			alterweight = Array{Int64}(undef,length(zll));
			maxweight = maximum(zll);
			for i = 1:length(alterweight)
				alterweight[i] = 1+maxweight-zll[i]
			end
			lowfilt = yafterx(alterweight,jrv)
			invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
			inversevec0 = Array{Int64}(undef,nl)
			inversevec0[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec0)
			jz = yafterx(jz,inversevec0)

			lowfilt = yafterx(ff2pv,jrv)
			invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
			inversevec1 = Array{Int64}(undef,nl)
			inversevec1[invertiblevec]=1:nl
			jrv = yafterx(jrv,inversevec1)
			jz = yafterx(jz,inversevec1)
			translatorvecb = yafterx(inversevec0,inversevec1)
			inversevec0 = [];inversevec1 = [];lowfilt = [];invertiblevec = [];
			#gc()
			(rt,ct,zt) = transposeLighter(jrv,jcp,jz,nll)
			colsum = ct.-1

			pmhist = zeros(Int64,m+1,m) #for sloth (apologies) we'll leave some unsed stuff in row m+1
			fpi = zeros(Int64,m+1,m)
			processfpi!(pmhist,fpi,jcp,jrv,ff2pv,m)

			#### reset ff2pv for next round
			ff2pvold = copy(ff2pv)
			ff2pv = Array{Int64}(undef,nl)
			ff2pv[:].= m+1

			oldclaw = Array{Int64}(undef,m)
		end

		for i = 1:m
			izfull[:].=0
			lrange = cran(jcp,i)
			izfull[jrv[lrange]] = jz[lrange]

			for j = (i+1):m
				dij = symmat[j,i]
				if dij == 0
					continue
				end
				if sd < maxsd-1
					process_sd_lt_maxsd!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
						npsupp::BitArray{1})
				elseif sd == maxsd-1
					process_sd_onelt_maxsd_1!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
						npsupp::BitArray{1})
				else
					for l = 1:(i-1)
						oldclaw[l] = minimum(symmat[l,[i,j]])
					end
					process_maxsd_one2i!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxsd_i2i!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxsd_i2j!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxsd_j2j!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})

					process_maxsd_j2end!(
						i::Int64,j::Int64,dij::Int64,stepsize::Int64,
						facecount::Array{Int64,1},numpairs::Array{Int64,1},
						jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
						r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
						oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
						rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
						izfull::Array{Int64,1},ff2pv::Array{Int64,1},
						pmhist::Array{Int64,2},fpi::Array{Int64,2},
						npsupp::BitArray{1})
				end
			end
			# update the column pattern and the total number of nonzeros
			# encountered per codim2 face
			c[i+1] = facecount[1]+1
			if sd == maxsd
				colsum[jrv[cran(jcp,i)]].+=1
			end
		end
		delrange = c[end]:length(r)
		deleteat!(r,delrange)
		deleteat!(z,delrange)
		deleteat!(pflist,(numpairs[1]+1):nl)
		if sd == maxsd
			r = translatorvecb[r]
		end
		firstv[sd] = c
		farfaces[sd] = r
		prepairs[sd] = pflist
		grain[sd] = z
		if isempty(farfaces[sd])
			for nextcard = (sd+1):maxsd
				firstv[nextcard] = [1;1]
				farfaces[nextcard] = Array{Int64}(undef,0)
				prepairs[nextcard] = Array{Int64}(undef,0)
				grain[nextcard] = Array{Int64}(undef,0)
			end
			if verbose
				println("no simplices of cardinality $(sd) or higher")
			end
			break
		end
	end
	if verbose
		println("collecting garbage")
		println(["number of edges" length(farfaces[2])])
	end
	#gc()
	if dictionaryoutput == true
		D = Dict{String,Any}(
			"farfaces" => farfaces,
			"firstv" => firstv,
			"grain" => grain,
			"prepairs" => prepairs,
			"symmat" => symmat,
			"nvl2ovl"=> vperm)
		return D
	else
		return farfaces,firstv,grain,prepairs,symmat,vperm
	end
end

function process_sd_lt_maxsd!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
	npsupp::BitArray{1})
	for k = cran(jcp,j)
		kk = jrv[k]
		farfilt = jz[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
			if claw >= farfilt && npsupp[k]
				pairupdate!(k,facecount,pflist,numpairs,npsupp,1)
			end
		end
	end
end

function process_sd_onelt_maxsd_1!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},pmhist::Array{Int64,2},
	npsupp::BitArray{1})
	for k = cran(jcp,j)
		kk = jrv[k]
		farfilt = jz[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
			if npsupp[k] && (claw >= farfilt)
				pairupdatedeluxe!(k,i,j,numpairs,facecount,pflist,ff2pv,npsupp,pmhist)
			end
		end
	end
end

function process_maxsd_one2i!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for l = 1:(i-1)
		if fpi[l,j]<fpi[l+1,j]
			ocl = oldclaw[l]
			if ocl < dij
				process_maxsd_one2i_subroutine!(
					i::Int64,j::Int64,dij::Int64,stepsize::Int64,
					facecount::Array{Int64,1},numpairs::Array{Int64,1},
					jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
					r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
					oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
					rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
					izfull::Array{Int64,1},ff2pv::Array{Int64,1},
					pmhist::Array{Int64,2},fpi::Array{Int64,2},
					npsupp::BitArray{1},
					l::Int64,ocl::Int64)
			end
		end
	end
end

function process_maxsd_one2i_subroutine!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1},
	l::Int64,ocl::Int64)
	for k = fpi[l,j]:(fpi[l+1,j]-1)
		kk = jrv[k]	## may have to reindex this
		farfilt = jz[k]
		if zll[kk] <= ocl
			break
		elseif oldclaw[l] < min(farfilt,izfull[kk])
			claw = min(izfull[kk],dij)
			if claw >= farfilt
				if npsupp[k]
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
					pairupdate!(k,facecount,pflist,numpairs,npsupp,3)
					ff2pv[k] = i
				elseif oldclaw[ff2pv[k]]>=farfilt
					continue
				elseif saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			elseif (claw>0) && saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,claw,stepsize)
			end
		end
	end
end

function process_maxsd_i2i!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[i,j]:(fpi[i+1,j]-1)
		kk = jrv[k]
		farfilt = jz[k]
		if dij >= farfilt && npsupp[k]
			faceupdate!(facecount,r,z,k,farfilt,stepsize)
			pairupdate!(k,facecount,pflist,numpairs,npsupp,4)
			ff2pv[k] = i
		else
			farfilt = min(dij,farfilt)
			if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
			end
		end
	end
end

function process_maxsd_i2j!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[i+1,j]:(fpi[j,j]-1)
		kk = jrv[k]
		if izfull[kk]>0
			farfilt = jz[k]
			claw = min(izfull[kk],dij)
			if claw >= farfilt && npsupp[k]
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
				pairupdate!(k,facecount,pflist,numpairs,npsupp,5)
				ff2pv[k] = i
			else
				farfilt = min(claw,farfilt)
				if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			end
		end
	end
end

function process_maxsd_j2j!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[j,j]:(fpi[j+1,j]-1)
		kk = jrv[k]
		if izfull[kk]>0
			claw = min(izfull[kk],dij)
			if saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
				faceupdate!(facecount,r,z,k,claw,stepsize)
			end
		end
	end
end

function process_maxsd_j2end!(
	i::Int64,j::Int64,dij::Int64,stepsize::Int64,
	facecount::Array{Int64,1},numpairs::Array{Int64,1},
	jrv::Array{Int64,1},jcp::Array{Int64,1},jz::Array{Int64,1},
	r::Array{Int64,1},z::Array{Int64,1},pflist::Array{Int64,1},
	oldclaw::Array{Int64,1},zll::Array{Int64,1},colsum::Array{Int64,1},
	rt::Array{Int64,1},ct::Array{Int64,1},zt::Array{Int64,1},
	izfull::Array{Int64,1},ff2pv::Array{Int64,1},
	pmhist::Array{Int64,2},fpi::Array{Int64,2},
	npsupp::BitArray{1})
	for k = fpi[j+1,j]:(jcp[j+1]-1)
		kk = jrv[k]
		if izfull[kk]>0
			farfilt = jz[k]
			claw = min(izfull[kk],dij)
			if claw >= farfilt && npsupp[k]
				faceupdate!(facecount,r,z,k,farfilt,stepsize)
				pairupdate!(k,facecount,pflist,numpairs,npsupp,6)
				ff2pv[k] = i
			else
				farfilt = min(claw,farfilt)
				if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
					faceupdate!(facecount,r,z,k,farfilt,stepsize)
				end
			end
		end
	end
end

function pairupdate!(k::Int64,facecount::Array{Int64,1},pflist::Array{Int64,1},numpairs::Array{Int64,1},npsupp::BitArray{1},iterateNumber)
	numpairs[1]+=1
	pflist[numpairs[1]] = facecount[1]
	npsupp[k]=false
end

function pairupdatedeluxe!(k::Int64,i::Int64,j::Int64,numpairs::Array{Int64,1},facecount::Array{Int64,1},pflist::Array{Int64,1},ff2pv::Array{Int64,1},npsupp::BitArray{1},pmhist::Array{Int64,2})
	numpairs[1]+=1
	pmhist[i,j]+=1
	npsupp[k]=false
	pflist[numpairs[1]]=facecount[1]
	ff2pv[k] = i
end

function faceupdate!(facecount::Array{Int64,1},r::Array{Int64,1},z::Array{Int64,1},k::Int64,farfilt::Int64,stepsize::Int64)
	facecount[1]+=1
	if facecount[1]>length(r)
		append!(r,Array{Int64}(undef,stepsize))
		append!(z,Array{Int64}(undef,stepsize))
	end
	r[facecount].= k
	z[facecount].= farfilt
end

function faceupdatedeluxe!(facecount::Array{Int64,1},r::Array{Int64,1},z::Array{Int64,1},k::Int64,farfilt::Int64,stepsize::Int64,s::Array{Int64,1},i::Int64)
	facecount[1]+=1
	if facecount[1]>length(r)
		append!(r,Array{Int64}(undef,stepsize))
		append!(z,Array{Int64}(undef,stepsize))
		append!(s,Array{Int64}(undef,stepsize))
	end
	r[facecount].= k
	z[facecount].= farfilt
	s[facecount].= i
end

function saveface(ct::Array{Int64,1},kk::Int64,colsum::Array{Int64,1},farfilt::Int64,oldclaw::Array{Int64,1},rt::Array{Int64,1},zt::Array{Int64,1})
	keep = true
	for l = ct[kk]:colsum[kk]
		if  zt[l]>= farfilt && oldclaw[rt[l]]>=farfilt
			keep = false
			break
		end
	end
	return keep
end


function processfpi!(pmhist::Array{Int64,2},fpi::Array{Int64,2},jcp::Array{Int64,1},jrv::Array{Int64,1},ff2pv::Array{Int64,1},m::Integer)
	for p = 1:m
		for q = jcp[p]:(jcp[p+1]-1)
			pmhist[ff2pv[jrv[q]],p]+=1
		end
	end
	for p = 1:m
		fpi[1,p] = jcp[p]
		for q = 1:m
			fpi[q+1,p] = fpi[q,p]+pmhist[q,p]
		end
	end
end

function generate2faces(symmat)
	m = size(symmat,1)
	if issparse(symmat)
		return symmat
	else
		L = 0
		for i = 1:m
			for j = (i+1):m
				if symmat[j,i]>0
					L+=1
				end
			end
		end
		rowval = Array{Int64}(undef,L)
		nzval = Array{Int64}(undef,L)
		colptr = Array{Int64}(undef,m+1)
		marker = 0
		colptr[1] = 1
		for i = 1:m
			colptr[i+1]=colptr[i]
			for j = (i+1):m
				if symmat[j,i]>0
					colptr[i+1]+=1
					rowval[colptr[i+1]-1] = j
					nzval[colptr[i+1]-1] = symmat[j,i]
				end
			end
		end
	end
	return rowval,colptr,nzval
end

function generate3faces!(
	farfaces_cell,
	firstv_cell,
	grain_cell,
	prepairs_cell,
	m,
	symmat;
	verbose = false)

	grain::Array{Int64,1} = grain_cell[2]
	farfaces::Array{Int64,1} = farfaces_cell[2]
	firstv::Array{Int64,1} = firstv_cell[2]

	numverts = length(firstv)-1
	numedges = length(farfaces)
	stepsize = size(symmat,1)^2
	facecount= [0]
	numpairs = 0

	closefaces = Array{Int64}(undef,numedges)
	for i = 1:m
		closefaces[cran(firstv,i)].=i
	end
	iso = integersinsameorder(farfaces)
	closefaces_higsorted = 	Array{Int64}(undef,numedges)
	grain_higsorted = 	Array{Int64}(undef,numedges)
	closefaces_higsorted[iso] = closefaces
	grain_higsorted[iso] = grain

	firstv_hs = zeros(Int64,m+1)
	for face in farfaces
		firstv_hs[face+1]+=1
	end
	firstv_hs[1] = 1
	for i = 2:(m+1)
		firstv_hs[i] = firstv_hs[i-1]+firstv_hs[i]
	end

	adist = Array{Int64}(undef,m)
	idist = Array{Int64}(undef,m)
	r = Array{Int64}(undef,numedges)
	z = Array{Int64}(undef,numedges)
	s = Array{Int64}(undef,numedges)

	clawvec = Array{Int64}(undef,m)
	ncheckedges = trues(numedges)

	for a = 1:m
		adist[:].=0
		adist[crows(firstv,farfaces,a)] = crows(firstv,grain,a)
		for ip = cran(firstv,a)
			i = farfaces[ip]
			dai = grain[ip]
			idist[:].=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
			for jp = cran(firstv,i)
				if ncheckedges[jp]
					j = farfaces[jp]
					dij = grain[jp]
					if dij <= dai && dij <= adist[j] # note this condition bakes in the req. that j be adjacent to a
						numpairs+=1
						ncheckedges[jp] = false
						clawvec[1:i] .= 0
						for lp = cran(firstv_hs,j)
							l = closefaces_higsorted[lp]
							if l >= i
								break
							elseif idist[l]!=0
								clawvec[l] = min(idist[l],grain_higsorted[lp])
							end
						end
						for kp = cran(firstv,j)
							k = farfaces[kp]
							djk = grain[kp]
							dak = adist[k]
							dik = idist[k]
							if dak < dij && dak<djk && dak<dik	# this bakes in req. that dik>0
								dijk = min(dij,dik,djk)
								keepface = true
								for bp = cran(firstv_hs,k)
									b = closefaces_higsorted[bp]
									if b >= i
										break
									elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
										keepface = false
										break
									end
								end
								if keepface
									faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
								end
							end
						end
					end
				end
			end
		end
	end
	holdi = 0
	for edge = findall(ncheckedges)
		i = closefaces[edge]
		j = farfaces[edge]
		dij = grain[edge]
		if i != holdi
			idist[:].=0
			idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
			idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
			holdi = i
		end
		clawvec[1:i] .= 0
		for lp = cran(firstv_hs,j)
			l = closefaces_higsorted[lp]
			if l >= i
				break
			elseif idist[l]!=0
				clawvec[l] = min(idist[l],grain_higsorted[lp])
			end
		end
		#### a facsimile of above
		for kp = cran(firstv,j)
			k = farfaces[kp]
			dik = idist[k]
			if dik==0
				continue
			end
			djk = grain[kp]
			dijk = min(dij,dik,djk)
			keepface = true
			for bp = cran(firstv_hs,k)
				b = closefaces_higsorted[bp]
				if b >= i
					break
				elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
					keepface = false
					break
				end
			end
			if keepface
				faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
			end
		end
		####
	end

	num3faces = facecount[1]
	holderlengths = length(r)
	deletionrange = (num3faces+1):holderlengths
	deleteat!(r,deletionrange)
	deleteat!(z,deletionrange)
	deleteat!(s,deletionrange)

	iso = integersinsameorder(s)
	r[iso]=r
	z[iso]=z
	fv3 = zeros(Int64,numverts+1)
	fv3[1] = 1
	for face in s
		fv3[face+1]+=1
	end
	for i = 2:(numverts+1)
		fv3[i] = fv3[i-1]+fv3[i]
	end

	pairmarker = 0
	npes = trues(numedges)
	prepairs = Array{Int64}(undef,numedges)
	for i = 1:num3faces
		edge = r[i]
		if npes[edge] && z[i] == grain[edge]
			npes[edge]=false
			pairmarker+=1
			prepairs[pairmarker]=i
		end
	end
	deleteat!(prepairs,(pairmarker+1):numedges)

	farfaces_cell[3]=r
	firstv_cell[3]=fv3
	grain_cell[3]=z
	prepairs_cell[3]=prepairs

	buffer1 = Array{Array{Int64,1},1}(undef,1)
	buffer2 = Array{Array{Int64,1},1}(undef,1)
	buffer1[1] = Array{Int64}(undef,0)
	buffer2[1] = ones(Int64,numverts+1)
	append!(farfaces_cell,buffer1)
	append!(grain_cell,buffer1)
	append!(prepairs_cell,buffer1)
	append!(firstv_cell,buffer1)

	return r,fv3,z,prepairs,numpairs
end


##########################################################################################

####	MAIN

##########################################################################################

function inputVmodel2defaultgeneraformat(s,model)
	if typeof(s) == String
		if in(model,["vr","pc"])
			entryformat = "textfile"
		elseif in(model,["complex"])
			entryformat = "sp"
		end
	else
		entryformat = "n/a"
	end
	return entryformat
end

"""

    eirene(X[, keyword arugemts])

Computes the persistent homology of a filtered complex.

"""
function eirene(
	s;
	model		= "vr",
	maxdim 		= 1,
	minrad		= -Inf,
	maxrad		= Inf,
	numrad		= Inf,
	nodrad 		= [],
	filfun	 	= "n/a",
	fastop		= true,
	vscale		= "default",
	record		= "cyclerep",
	entryformat  = inputVmodel2defaultgeneraformat(s,model),
	pointlabels	= [],
	verbose		= false)

	if in(model,["vr","pc"])
		maxsd = 		maxdim+2
		D = persistf2vr(
			s,
			maxsd;
			model 		= model,
			minrad 		= minrad,
			maxrad 		= maxrad,
			numrad 		= numrad,
			nodrad 		= nodrad,
			fastop 		= fastop,
			filfun 		= filfun,
			record 		= record,
			entryformat= entryformat,
			pointlabels = pointlabels,
			verbose 	= verbose)

		return 	D
	elseif model == "complex"
		D   = 	persistf2complex(
				s;
				maxdim=maxdim,
				entryformat=entryformat,
				record = record,
				verbose=false)
		return 	D
	else
		println()
		println("Error: the only valid values for keyword <model> are \"vr\", \"pc\", and \"complex\".")
		println("user input:")
		println(model)
	end
end

function genera2autoformat(rv,dp,dv,ev)
	if 	typeof(rv) 	== 	Array{Array{Int64}}
		return "segmented complex"
	elseif 	!isempty(dp)
		return "dp"
	elseif 	!isempty(dv)
		return "dv"
	elseif 	!isempty(ev)
		return "ev"
	end
end

#=
NB: the default value for maxdim has not been tested, and may cause errors;
the -3 accounts for (1) julia uses 1-indexed arrays, (2) to calculate
homology in dimension p, one must inspect the (p+1)-dimensional boundary
operator, (3) this operator should be givent the same treatment as those
that precede it ... generally this assumes that the next one up is at least
defined, even if it is trivial.
=#
function eirene(   	;
	rv				= 	zeros(Int64,0),
	cp 				= 	zeros(Int64,0),
	fv				= 	zeros(Int64,0),
	dp 				= 	zeros(Int64,0),
	dv 				= 	zeros(Int64,0),
	ev 				= 	zeros(Int64,0),
	model 			= 	"complex",
	entryformat 	= 	genera2autoformat(rv,dp,dv,ev),
	filfun 			= 	"n/a",
	maxdim 			= 	[],
	record			=	"cyclerep",
	pointlabels		=	[],
	verbose			=	false)

	D = persistf2complex(
		rv			= 	rv,
		cp			= 	cp,
		fv 			= 	fv,
		dp 			= 	dp,
		dv 			= 	dv,
		ev		= 	ev,
		maxdim 		= 	maxdim,
		record 		= 	record,
		verbose 	= 	verbose)

	return D
end


function checkdv(rv_ag,cp_ag,dv)
	for p 				=	1:length(cp_ag)-1
		if 	any(dv[crows(cp_ag,rv_ag,p)].!=(dv[p]-1))
			println("error: please check dimension values (dv)")
			return [rv_ag,cp_ag, dv]
		end
	end
	return []
end



function eulervector2dimensionpattern(ev)
	l 				= 	length(ev)
	dp 				= 	zeros(Int64,l+1)
	dp[1]			= 	1
	for p 			= 	1:l
		dp[p+1] 	= 	dp[p]+ev[p]
	end
	return 				dp
end


function dimensionvalues2dimensionpattern(dv)
	m 				= 	maximum(dv)
	dp 				= 	zeros(Int64,m+2)
	dp[1] 			= 	1
	currentsd 		= 	1
	for p 			= 	1:length(dv)
		while dv[p] > 	currentsd-1
			currentsd			+= 	1
			dp[currentsd]		= 	p
		end
	end
	dp[end] 		= 	length(dv)+1
	return dp
end

function offdiagmin(d::Array{Tv}) where Tv
	if 	size(d,1) != size(d,2)
		println()
		println("error: d should be square")
		return
	end
	v 		= 	zeros(Tv,size(d,2))
	for 	p 	= 	1:size(d,2)
		val1 	= 	empteval(minimum,d[1:p-1,p],Inf)
		val2 	= 	empteval(minimum,d[p+1:end,p],Inf)
		v[p] 	= 	min(val1,val2)
	end
	return 	v
end

function ceil2grid(M;origin=0,stepsize=1,numsteps=Inf)
	if 	stepsize 	<  	0
		println()
		println("error in function <roundentries>: stepsize must be positive")
		return
	end
	if 	numsteps 	<  	1
		println()
		println("error in function <roundentries>: numsteps must be positive")
		return
	end

	N 				= 	copy(M)
	N 				= 	Array{Float64}(N) # conversion
	N  				= 	(N .- origin)./stepsize
	N 				= 	ceil.(N)
	N 				= 	N.*stepsize.+origin

	N[N.<origin]	.= 	-Inf
	if 	numsteps 	< 	Inf
		maxval 		= 	origin+numsteps*stepsize
		N[N.>maxval].= 	Inf
	end
	return 			N
end

###### Adding Wasserstein distances between persistence diagrams ##########

include("wasserstein_distances.jl")
#print("included wasserstein")
end # module
