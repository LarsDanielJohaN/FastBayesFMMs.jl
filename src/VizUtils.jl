
using Statistics,  Base.Threads  , Plots, NaNStatistics, Statistics ,  KernelDensity, Base.Threads


#=
postf, compute bounds for (1-alpha)% credible bands per points in design grid and mean function. 

F: Matrix{Float64} of dimensions Tn x Nf where Tn is the number of points in the design grid and Nf is the number of functions. 
alpha: Float64 or Nothing, significance level, defaults to 0.05. 
=#

function postf(F, alpha::Union{Float64, Nothing} = 0.05)
    @assert (alpha >= 0.0) & (alpha <= 1.0) "Fatal error on postf, alpha should satisfy 0.0 <= alpha <= 1.0" 

    Tn, Nf = size(F) #Obtains dimensions of F
    q_alpha_half= zeros(Tn) #Declares vector of dimension Tn to store empiricial quantiles at alpha/2 
    q_one_m_alpha_half = zeros(Tn) #Declares vector of dimension Tn to store empirical quantiles at 1 - alpha/2
    mu = mean(F, dims = 2) #Computes point wise means. 

    Threads.@threads for tn in 1:Tn #Computes quantiles per design point in parallel. 
        q_alpha_half[tn] = quantile( F[tn, :], alpha/2)
        q_one_m_alpha_half[tn] = quantile( F[tn, :], 1- alpha/2)
    end

    return Dict("mu"=>mu, "q_alpha_half"=>q_alpha_half, "q_one_m_alpha_half"=>q_one_m_alpha_half)

end



function plotf_mean_q(F, T_grid::Vector{Float64} = range(0, 1, length = size(F)[1]) |> collect,  alpha::Float64 = 0.05, alpha_viz::Float64 = 0.5, add_area::Bool = true)
    @assert length(T_grid) == size(F)[1] "Fatal error on plotf! T_grid should be of size Tn. "

    dict_postf = postf(F, alpha) #Obtains pointwise quantile values and mean function. 
    pq = plot(T_grid, dict_postf["q_one_m_alpha_half"], fillrange = dict_postf["q_alpha_half"], fillalpha = alpha_viz, color = :red, label = false, linealpha = 0) #Adds the filling between confidence bands.
    pq =plot!(pq, T_grid, [dict_postf["mu"], dict_postf["q_alpha_half"], dict_postf["q_one_m_alpha_half"]], label = ["mu" "q: $(alpha/2)" "q: $(1-alpha/2)"], color = [:red :blue :blue]) #Adds confidence bands. 
    return pq 


end
#=
plotf returns a plot object with a given collection of functions in it. 
F: Matrix{Float64} of dimensions Tn x Nf where Tn is the number of design points and Nf the total number of functions to plot. 
T_grid: Union{Vector{Float64}, Nothing}, vector of length Tn with design points. Defaults to an equally spaced grid of points in the unit interval of size Tn. 
add_postf: Union{Bool, Nothing}, boolean value which tells whether to add confidence bands and mean value, defaults to true. If add_postf = true, the functions in F 
will be ploted with a lower intensity than the mean and confidence bands. 
alpha: Union{Float64, Nothing} significance level for the confidence bands to use in case add_postf = true, defaults to 0.05. 
alpha_viz: Union{Float64, Nothing}, intensity of curve values. e.g. the alpha value specified to plot. 
=#

function plotf(F::Matrix{Float64}, T_grid::Union{Vector{Float64}, Nothing} = range(0, 1, length = size(F)[1]) |> collect, add_postf::Union{Bool, Nothing} = true, alpha::Union{Float64, Nothing} = 0.05, alpha_viz::Union{Float64, Nothing} = 0.2)
    @assert length(T_grid) == size(F)[1] "Fatal error on plotf! T_grid should be of size Tn. "
    Tn, Nf = size(F)
    p = plot(T_grid, F, label = false, alpha = alpha_viz) #Plots functions on p
    if add_postf #If add_postf = true, then obtain confidence bands, mean function and plot them. 
        dict_postf = postf(F, alpha)
        p =plot!(p, T_grid, [dict_postf["mu"], dict_postf["q_alpha_half"], dict_postf["q_one_m_alpha_half"]], label = ["mu" "q: $(alpha/2)" "q: $(1-alpha/2)"], color = [:red :blue :blue])
    end

    return p

end


#T_grid: Union{Vector{Float64}, Nothing}, vector of length Tn with design points. Defaults to an equally spaced grid of points in the unit interval of size Tn. 
#flfosr_output: Dict{String,Any}. Dictionary designed to take the output of flfosr as input. 
function compute_plot_residuals(flfosr_output::Dict{String,<:Any}; T_grid:: Union{Vector{Float64}, Nothing} = range(0, 1, length = size(Y)[1]) |> collect, alpha::Union{Float64, Nothing} = 0.3, title::Union{String, Nothing} = "")
    @assert keys(flfosr_output) == keys(Dict("X"=>1, "B"=>1,  "w_post"=>1, "ga_post"=>1, "alpha_post"=>1, "alpha_postf"=>1, "sig_eps_post"=>1, "sig_alpha_post"=>1, "sig_gamma_post"=>1, "sig_omega_post"=>1, "Y_hat" =>1, "Y_residuals"=>1 )) "Fatal error on compute_plot_residuals! Dictionary names for flfosr_output must match\nthose of an output dictionary from the flfosr method."
    println("Warning! Assertions only include input types and \nmatching names of a provided dictionary and those which are given by a flfosr output.\nContents run at personal risk.") 

    residuals_v = flfosr_output["Y_residuals"] #Gets residuals from all iterations.
    mean_T = mean(residuals_v, dims=1) #Calculates residual estimates throught iterations.
    aux_T = mean_T[1, :, :] #Gets corresponding mean estimate matrix. 
    tot_plot = scatter(T_grid, aux_T, alpha = alpha, color = "gray" , legend = false) #Creates scatter plot of domain vs residuals. 
    hline!(tot_plot, [0], color=:red, linestyle=:dash, lw = 1.5, label="y = 0") #Adds a horizontal line at y = 0. 

    return tot_plot 
end

#Mean residual kern_estimate. 
#Calculates the histogram of mean residuals. 
#i.e. of the MCMC residuals, calculate the mean and obtain a kernel estimate.


function mean_residual_kern_estimate(flfosr_output::Dict{String,<:Any}; no_points_grid::Union{Int64, Nothing} = 2000, color::Union{String, Nothing} = "black", alpha::Union{Float64, Nothing} = 0.3, title::Union{String, Nothing} = "", legend_title::Union{String, Nothing}  = "")
    @assert keys(flfosr_output) == keys(Dict("X"=>1, "B"=>1,  "w_post"=>1, "ga_post"=>1, "alpha_post"=>1, "alpha_postf"=>1, "sig_eps_post"=>1, "sig_alpha_post"=>1, "sig_gamma_post"=>1, "sig_omega_post"=>1, "Y_hat" =>1, "Y_residuals"=>1 )) "Fatal error on compute_plot_residuals! Dictionary names for flfosr_output must match\nthose of an output dictionary from the flfosr method."
    println("Warning! Assertions only include input types and \nmatching names of a provided dictionary with those which are given by a flfosr output.\nContents run at personal risk.")

    residuals_v = flfosr_output["Y_residuals"] #Gets residuals from all iterations. 
    mean_T = mean(residuals_v, dims=1) #Calculates residuals estimates throught iterations. 
    aux_T = mean_T[1, :, :] #Gets corresponding mean estimate matrix. 
    kde_mean_res = kde(vec(aux_T)) #Creates kernel estimate with all of residual estimates. 
    #Plots kernel estimate. 
    tot_plot = plot(range(min_val, max_val, no_points_grid ) |> collect, kde_mean_res.x, kde_mean_res.density, color = color, alpha = alpha, title = title,  )
    hline(tot_plot, [0], color =:red, linestyle = :dash, lw = 1.5, label = "x = 0") #Adds a vertical line at x = 0, should be mode if assumptions are met. 
    return tot_plot
end


#Calculate the variance of a posterior estimand. 
#C: A N x M matrix where N is the length of each MCMC chain and M is the number of chains. 
function posterior_estimand_convergenceq(C :: AbstractMatrix{Float64})
    N,M = size(C) #Gets the dimensions of C
    @assert (N > 3) & (M >= 2) "Fatal error on posterior_estimand_variance!!  Either (steps per chain) N>=4 or (number of chains) M >= 2 isnt fullfiled.  "
    #println("Warning!! It is only asserted that (number of steps per chain) N >= 4 and (number of chains) M >= 2. The rest goes at personal risk. ")
    chain_means = mean(C, dims = 1) #Calculate M means for each of the chains. 
    B = N*var(chain_means) #Calculates between sequence variance (B) term from BDA3.
    s_j_all = var(C, dims = 1) #Calculates M within chain empirical variances. 
    W = mean(s_j_all) #Calculates the mean within-sequence variances. 
    var_post = (W*(N-1) + B)/N #Calculates posterior variance of the estimand (an overestimate in terms of BDA3).
    R = sqrt(var_post/W) #Calculates estimation for potential scale reduction. 


    rho_pot = Vector{Float64}(undef, N-1) #Gets array for all potential rho values.
    Vt = 0.0 #assings place for variogram calculations. 

    for t in 1:(N-1)
        #mapreduce applies function (x,y) -> f(x,y) per pair elements A[i] and B[i] 
        Vt = mapreduce((x,y)-> abs2(x-y), +, view(C, (t+1):N, :), view(C, 1:(N-t), :) )/(M*(N-t)) #Gets variogram for lag t. 
        rho_pot[t] = 1 - Vt/(2*var_post) #Gets autocorrelation estimate of sequence at lag t. 
    end 

    use = Int64(floor((N-2)/2) ) 

    aux_use = view(rho_pot, 2:2:(2*use)) .+ view(rho_pot, 3:2:(2*use +1)) #Computes p_{2t} + p_{2t + 1}
    cut_off = findfirst(i -> (aux_use[i]<0) && isodd(i), eachindex(aux_use))

    rho_sum = if cut_off === nothing #If there isnt an odd T such that p_{T+1} + p_{T+2}<0, add everything except possibly the last element. 
        sum(aux_use) + rho_pot[1]
    elseif cut_off > 1 #If there is such a T and p_{2} + p_{3}<0 (i.e. there are other elements to add) add everything up to the cutoff
        sum(view(aux_use, 1:(cut_off - 1))) + rho_pot[1]
    else #In other case since p_{2} + p_{3} < 0, you only return the first element. 
        rho_pot[1]
    end
    n_eff = (M*N)/(1 + 2*rho_sum) #Computes number of effective samples. 

    return Dict("var_post"=>var_post, "B"=>B, "W"=>W, "R" => R, "n_eff"=>n_eff, "lag_correlations"=>rho_pot) #Returns a dictionary with the requried elements. 

end


function assess_coef_convergence(multiple_flfosr_output::Vector{Dict{String, Array{Float64}}} )
    @assert length(multiple_flfosr_output) > 1 "Fatal error on assess_coef_convergence! This functions requires at least two outputs to flfosr to assess convergence!!"
    ref_keys = keys(Dict("X"=>1, "B"=>1,  "w_post"=>1, "ga_post"=>1, "alpha_post"=>1, "alpha_postf"=>1, "sig_eps_post"=>1, "sig_alpha_post"=>1, "sig_gamma_post"=>1, "sig_omega_post"=>1, "Y_hat" =>1, "Y_residuals"=>1 )) 
    @assert all(d -> keys(d) == ref_keys, multiple_flfosr_output) "Fatal error on assess_coef_convergence! At least one of the dictionaries in multiple_flfosr_output has keys different from those produced by flfosr. "
    println("Warning! Assertions only include input types and \nmatching names of a provided dictionary with those which are given by a flfosr output.\nContents run at personal risk.")
    println("Important note: For large L and Tn and Niter MCMC iterations per function, there will be (L+1)*Tn on which to check convergence.\n Thus this function may take a while to run. ")

    Tn, Lp1, Niter = size(multiple_flfosr_output[1]["alpha_postf"]) #Gets the dimensions required. 
    Alpha_n_eff = Matrix{Float64}(undef, Tn, Lp1) #Gets matrix to store the number of effective samples. 
    Alpha_n_R = Matrix{Float64}(undef, Tn, Lp1) #Gets matrix to store the number of R values. 

    @views begin 
        Threads.@threads for l in 1:Lp1
            for tn in 1:Tn 
                aux = posterior_estimand_convergenceq(hcat((multiple_flfosr_output[i]["alpha_postf"][tn, l, :] for i in 1:length(multiple_flfosr_output))...))
                Alpha_n_eff[tn, l] = aux["n_eff"]
                Alpha_n_R[tn, l] = aux["R"]
            end
        end
    end

    return Dict("function_point_n_eff"=>Alpha_n_eff, "function_point_R"=>Alpha_n_R)
end 


