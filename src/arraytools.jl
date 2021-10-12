"""
$(TYPEDSIGNATURES)

Resurrect linspace despite https://github.com/JuliaLang/julia/pull/25896#issuecomment-363769368
"""
linspace(a,b,n)=collect(range(a,b,length=n))


"""
$(TYPEDSIGNATURES)

(Try to) create a subdivision of interval (a,b) stored in the 
returned array X such that 
  - `X[1]==a, X[end]==b`
  - `(X[2]-X[1])<=ha+tol*(b-a)`
  - `(X[end]-X[end-1])<=hb+tol*(b-a)`
  - There is a number q such that  `X[i+1]-X[i] == q*(X[i]-X[i-1])`
  - X is the array with the minimal possible number of points with the above property
  
Caveat: the algorithm behind this is  tested for many cases but unproven.

Returns an Array containing the points of the subdivision.
"""
function geomspace(a, b, ha, hb ; tol=1.0e-10, maxiterations=100)
    a=Float64(a)
    b=Float64(b)
    ha=Float64(ha)
    hb=Float64(hb)
    
    function _geomspace0(l,h0, hl, tol=1.0e-10)
        @assert (l>0.0)
        @assert (h0>0.0)
        @assert (hl>=h0)
        @assert((hl+h0)<l)
        
        #  We need to  adjust two things:
        
        # The sum of the geometric progression must
        # match the length of the interval, so lmismatch
        # should be zero:
        function lmismatch(q,k)
            return l - h0*(1-q^k)/(1-q)
        end
        
        # The claim from experimenral evidence (Wolfram) 
        # is that, if written as a polynomial,
        # it has two real zeros: one and the value searched for which
        # is slightly larger than one. All other zeros are on one circle.
        
        # The size of the last interval should be close to 
        # to hl, so hmismatch should be close to one and not larger than one
        function  hmismatch(q,k)
            return  h0*q^(k-1)/hl
        end

        # define initial number of intervals from
        # average of minmal and maximal h
        n=Int(ceil((2.0*l/(h0+hl))))
        

        if n==1
            n=2
        end
        
        # define initial q such that hmismatch is one.
        q=(hl/h0)^(1.0/(n-1.0))
        
        # Iteration until both mismatches are satisfactory
        # Outer loop runs until hmismatch is less than 1
        hmiss=10.0 # some initial value >1 just to run the loop at least once
        if abs(q-1.0)<tol
            hmiss=1.0
        end

        while  hmiss>1.0
            # increase number of intervals until
            # lmismatch becomes less than zero
            while  lmismatch(q,n)>0.0
                n+=1
            end

            # find initial interval for q containing
            # value with zero lmismatch 
            ns=0
            
            while lmismatch(q,n)<0 &&  ns<maxiterations
                q*=0.99
                ns+=1
            end
            
            xl=q
            xr=q/0.99
            @assert ns<maxiterations "Unable to determine geomspace data after $(maxiterations) iterations"
            
            # bisection to define q with zero lmismatch
            ns=0
            xm=0.5*(xl+xr)
            while abs(xr-xl)>tol && ns<maxiterations
                ns+=1
                mmm=lmismatch(xm,n)
                if abs(mmm)<tol
                    break
                elseif   lmismatch(xl,n)*mmm<0
                    xr=xm
                else
                    xl=xm
                end
                xm=0.5*(xl+xr)
            end
            # increase q slightly to increase probability
            # for last interval to be <=hl
            q=xm*(1.0+100.0*tol)
            @assert ns<maxiterations
            hmiss=hmismatch(q,n)
            if hmiss>1.0+tol 
                n=n+1
            end
        end
        #  printf("%d %g %g %g\n",n,q,lmismatch(q,n),hmismatch(q,n))

        X = Array{Float64,1}(undef,n+1)
        X[1]=0
        h=h0
        for i=1:n
            X[i+1]=X[i]+h
            h*=q
        end
        X[n+1]=l

        # sometimes, we get a zero sized last interval
        # so we just pop the last element.
        if X[n]>X[n+1]
            pop!(X)
            n=n-1
            @assert X[n+1]-X[n] > hl/2 "Interval turned out to be less than $(hl)/2"
        end
        return X
    end

    @assert (ha>0.0) "Start step size $(ha) should be positive"
    @assert (hb>0.0) "End step size $(hb) should be positive"
    @assert (a<b)    "Interval ends $(a), $(b) should be increasing"
    @assert ((ha+hb)<b-a) "Sum of step sizes $(ha)+$(hb) should not exceed interval size $(b)-$(a)"

    
    # Map things to [0,b-a]
    tol=tol*(b-a)
    if ha<hb-tol
        X=_geomspace0(b-a,ha,hb,tol)
        X.+=a
    elseif ha>hb+tol
        X=-reverse(_geomspace0(b-a,hb,ha,tol))
        X.+=b
    else
        n=Int(ceil((b-a)/ha))
        X=collect(range(a,b,length=n+1))
    end
#    @show X[2]-X[1],ha, X[end]-X[end-1],hb

    @assert (X[2]-X[1])<=ha+tol  "First interval turned out to be larger than $(ha)"
    @assert (X[end]-X[end-1])<=hb+tol  "Last interval turned out to be larger than $(hb)"
    @assert X[begin]≈a  "Range start differs from $(b)"
    @assert X[end]≈b "Range end differs from $(b)"
    
    return X
end



collect_or_assign(X)=X
collect_or_assign(X::AbstractRange)=collect(X)
is_monotone(X)=all(X[2:end]-X[1:end-1].>0)

"""
    c=glue(a,b)

Glue together two vectors `a` and b resulting in a vector c. They last element 
of `a` shall be equal (up to tol) to the first element of b.
The result fulfills `length(c)=length(a)+length(b)-1`
"""
function glue(_a::AbstractVector, _b::AbstractVector; tol=1.0e-10)
    a=collect_or_assign(_a)
    b=collect_or_assign(_b)

    is_monotone(a) || error("non-monotonous first argument of glue")
    is_monotone(b) || error("non-monotonous second argument of glue")
    Tv=promote_type(eltype(a),eltype(b))
    na=length(a)
    nb=length(b)
    
    d=b[1]-a[na-1]
    @assert(d>0)
    d=b[1]-a[na]
    @assert(d>-tol)
    @assert(d<tol)

    c=Vector{Tv}(undef,na+nb-1)
    ic=0
    for ia=1:na
        ic+=1
        c[ic]=a[ia]
    end
    for ib=2:nb
        ic+=1
        c[ic]=b[ib]
    end
    return c
end
