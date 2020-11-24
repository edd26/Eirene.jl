
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

