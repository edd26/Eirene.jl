
function ocfcheckfun3()
	n = 10
	m = 1000
	for q = 1:n
		M = rand(m,m)
		M = M+M'
		if q < n-5
			numrad = rand(10:90)
			upperlim = rand()*2
		elseif q == n-5
			numrad = 10
			upperlim = rand()*2
		elseif q == n-4
			numrad = 1
			upperlim = -1
		elseif q == n-3
			numrad = 1
			upperlim = rand()*2
		elseif q == n-2
			numrad = 1
			upperlim = 0
		elseif q == n-1
			numrad = 1
			upperlim = 3
		elseif q == n
			numrad = 1
			upperlim = -1
		end

		ocf1 = Array{Any,1}(8)
		ocf2 = Array{Any,1}(8)
		j = 0
		for a in [Inf, upperlim]
			for b in [Inf,numrad]
				for c in [true, false]
					j+=1
					print("RRR1[[j=[$(j)]]]")
					ocf1[j] = ordercanonicalform_3(M,maxrad=a,numrad=b,fastop=c)
				end
			end
		end
		print("ocf3-done-")

		if numrad > 1
			N = copy(M)
			for i = 1:m
				N[i,i] = Inf
			end
			minentry = minimum(N)
			for i = 1:m
				N[i,i] = -Inf
			end
			maxentry = maximum(N)

			Q = minmaxceil(copy(M),minentry,upperlim,numrad)
			N = minmaxceil(copy(M),minentry,maxentry,numrad)
			for i = 1:m
				N[i,i] = Inf
				Q[i,i] = Inf
			end
		end
		if numrad == 1
			N = copy(M)
			if upperlim == Inf
				for i = 1:m
					N[i,i] = -Inf
				end
				maxrad = maximum(N)
			else
				maxrad = upperlim
			end
			N[N.<=maxrad]=maxrad
			N[N.>maxrad]= maxrad+1
			for i = 1:m
				N[i,i] = Inf
			end
			Q = copy(N)
		end

		j = 0
		for a in [Inf, upperlim]
			for b in [Inf,numrad]
				if b == numrad
					if a == upperlim
						L = copy(Q)
					else
						L = copy(N)
					end
				else
					L = copy(M)
				end

				for c in [true, false]
					j+=1
					LL = copy(M)

					for i = 1:m
						LL[i,i] = Inf
					end
					publicmin = minimum(LL)

					for i = 1:m
						LL[i,i]=-Inf
					end
					if a == Inf
						publicmax = maximum(LL)
					else
						publicmax = a
					end

					if publicmax < publicmin
						ocf2[j] = zeros(Int64,m,m),Array{Float64}(undef,0)
					elseif b == 1 && c
						privatemax = minimum(maximum(LL,1))
						if publicmax >= privatemax
							ocf = zeros(Int64,m,m)
							index = findfirst(maximum(LL,1),privatemax)
							ocf[index,:]=1
							ocf[:,index]=1
							ocf[index,index]=0
							ocf2[j]= (ocf,[publicmax],M)
						else
							ocf = zeros(Int64,m,m)
							ocf[LL.<=publicmax]=1
							for i=1:m
								ocf[i,i]=0
							end
							ocf2[j] = (ocf,[publicmax],M)
						end
					elseif b == 1 && a == Inf
						ocf = ones(Int64,m,m)
						for i=1:m
							ocf[i,i]=0
						end
						ocf2[j] = (ocf,[maximum(LL)],M)
					else
						ocf2[j] = ordercanonicalform(L,maxrad=a,fastop=c)
					end
				end
			end
		end
		for k = 1:8
			if ocf1[k][1]!=ocf2[k][1] || sum(abs(ocf1[k][2]-ocf2[k][2])) > 0.00000001
				println("locus = [$(k)]")
				println("q = [$(q)]")
				println("numrad = $(numrad)")
				println("upperlim = $(upperlim)")
				return ocf1,ocf2,M,N
			end
		end
	end
	return "passedtest"
end
