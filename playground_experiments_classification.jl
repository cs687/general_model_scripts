### Preparations ###

# Load functions.
include("../Codebase/Functions_specific/trajectory_classification.jl")

# Load data (simulated and experimental).
exp_trajss_all = load_all_trajectories()
exp_trajss = [[filter(e -> length(e) == 572, etjs) for etjs in exp_trajs] for exp_trajs in exp_trajss_all]
sim_trajss = load_all_simulations()

# Classifies first set of experiments.

# Plot classification examples.
get_info = (t,s)->""
plot_classification_examples(exp_trajss[1]; info_f=get_info)
plot_classification_examples(exp_trajss[2]; info_f=get_info)
plot_classification_examples(exp_trajss[3]; info_f=get_info)
plot_classification_examples(exp_trajss[4]; info_f=get_info)
plot_classification_examples(exp_trajss[5]; info_f=get_info)

sim_traj_examples = rand(sim_trajss)
plot_classification_examples(sim_traj_examples; info_f=get_info)
 
set_data = get_set_data(sim_traj_examples)
traj = sim_traj_examples[60][60]
plot(traj)
class = classify_trajectory_alt(traj, set_data)

@time plot_classification_examples(sim_traj_examples; info_f=get_info)
### Alternative Classifier ###
sd1 = get_set_data(exp_trajss[1])
sd1 = get_set_data(exp_trajss[2])
sd1 = get_set_data(exp_trajss[3])
sd1 = get_set_data(exp_trajss[4])
sd1 = get_set_data(exp_trajss[5])

set_data

dens = kde(vcat(vcat(exp_trajss[1]...)...))
pk_idxs, pk_vals = findmaxima(dens.density)
plot(dens.x, dens.density)


set_data = sd1
trajectory = exp_trajss[1][1][5]

plot(trajectory)
plot!([1, length(trajectory)], [set_data[:min_peak], set_data[:min_peak]]; color=:black, lw=2)
plot!([1, length(trajectory)], [set_data[:max_peak], set_data[:max_peak]]; color=:black, lw=2)

get_transition_types(trajectory, set_data)
classify_trajectory_alt(trajectory, set_data)

maximum(trajectory)
set_data[:min_peak]

density(vcat(vcat(exp_trajss[1]...)...))
plot!([dens.x[i], dens.x[i]],[0.0, 0.00100])
plot!([dens.x[i], dens.x[i]],[0.0, 0.00100])
plot!([dens.x[i], dens.x[i]],[0.0, 0.00100])
plot!([dens.x[i], dens.x[i]],[0.0, 0.00100])

using KernelDensity
dens = kde(vcat(vcat(exp_trajss[1]...)...))
dens |> print_fields

alls = findall((dens.density[2:end-1] .> dens.density[3:end]) .&& (dens.density[2:end-1] .> dens.density[1:end-2]))

# Activation time.
act_times = []
for (idx,val) in enumerate(traj[1:end-1])
    (val > t1) && continue
    act_time = findfirst(traj[idx+1:end] .> t2)
    isnothing(act_time) || push!(act_times, act_time)
end
act_times



function transition_times(traj, t1, t2)

end





### Temporary Analysis ###
# Helper function.
using LsqFit
begin
    function find_freq(wpd)
        idx1 = findfirst(wpd.power[2:end] .> wpd.power[1:end-1])
        isnothing(idx1) && (return wpd.freq[1])
        idx2 = idx1 + argmax(wpd.power[idx1:end]) - 1
        isnothing(idx2) && (return wpd.freq[1])
        wpd.freq[idx2]
    end
    function find_period(trajectory)
        wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 100)
        return round(1/(12*find_freq(wpd)), sigdigits=3)
    end
    # function plot_wpd(wpd)
    #     plot(wpd.freq,wpd.power)
    # end
    function plot_wpd(traj::Vector; n=200)
        wpd = welch_pgram(traj ./ mean(traj) .-1, n)
        plot(wpd.freq,wpd.power)
    end
    function plot_wpd_fit(traj::Vector; n=200)
        wpd = welch_pgram(traj ./ mean(traj) .-1, n)
        plot(wpd.freq,wpd.power)
        @. model(x, p) = p[1]*exp(-x*p[2])
        @. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
        fit = curve_fit(model, wpd.freq, wpd.power, [0.5, 0.5])
        plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.5))
    end
    function plot_wpd_fit_pink(traj::Vector; n=200)
        wpd = welch_pgram(traj ./ mean(traj) .-1, n)
        plot(wpd.freq,wpd.power)
        @. model(x, p) = p[1]*1/(x^p[2])
        fit = curve_fit(model, wpd.freq, wpd.power, [0.5, 0.5])
        plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.5))
    end
    function analyse_wpd_fit(traj::Vector; n=200)
        wpd = welch_pgram(traj ./ mean(traj) .-1, n)
        plot(wpd.freq,wpd.power; la=0.8,lw=7, label="Periodogram")
        @. model(x, p) = p[1]*exp(-x*p[2])
        @. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
        fit = curve_fit(model, wpd.freq, wpd.power, [0.5, 0.5])
        plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.1); la=0.8,lw=7, label="Exponential fitted to periodogram")

        proms = [calc_prominensce(wpd.power[i], model(wpd.freq[i], fit.param); pen=maximum(wpd.power)/20.0) for i = 1:length(wpd.freq)]
        proms = maximum(wpd.power) * proms ./ maximum(proms)
        plot!(wpd.freq, proms; la=0.8,lw=7, label="Periodogram prominensce")

        plot_freq_line!(1, maximum(wpd.power))
        plot_freq_line!(2, maximum(wpd.power))
        plot_freq_line!(4, maximum(wpd.power))
        plot!(framestyle=:box, xguide="Frequency",guidefontsize=46)
    end
    function calc_prominensce(v1, v2; pen=1.0)
        #return max(0.0, (v1 - v2))
        (v1 - pen < v2) && return 0.0
        return (v1 - pen) / v2
    end
    function plot_traj_grid(traj)
        plot([i/6 for i in 1:length(traj)], traj, xticks=[i/6 for i in 1:length(traj)][6:12:end],grid=true,gridalpha=1.0)
    end
    function plot_freq_line!(period, max)
        freq = 1/(12*period)
        plot!([freq,freq],[0.0,max], lw=3, color=:grey)
    end
    function get_peak_periods(traj; n=200)
        wpd = welch_pgram(traj ./ mean(traj) .-1, n)
        @. model(x, p) = p[1]*exp(-x*p[2])
        @. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
        fit = curve_fit(model, wpd.freq, wpd.power, [0.5, 0.5])
        proms = [calc_prominensce(wpd.power[i], model(wpd.freq[i], fit.param); pen=maximum(wpd.power)/20.0) for i = 1:length(wpd.freq)]
        idxs = find_peak_indexes(proms)
        periods = [1/(12*wpd.freq[idx]) for idx in idxs]
        return periods
    end
end

function find_peak_indexes(vals)
    for idx = 2:length(vals)-1
        (vals[idx] < vals[idx-1]) && (vals[idx] = 0.0)
        (vals[idx] < vals[idx+1]) && (vals[idx] = 0.0)
    end
    idxs = []
    for i = 1:1
        all(vals .== 0.0) && break
        push!(idxs, argmax(vals))
        vals[idxs[end]] = 0.0
    end
    filter!(idx -> 1<idx<length(vals), idxs)
    return idxs
end

# Plot period distributions.
plot()
osc_trajs = get_trajectory_class(exp_trajss[5], :oscillation)
periods = find_period.(osc_trajs)
osc_p_peaks = vcat(get_peak_periods.(osc_trajs)...)
density!(osc_p_peaks; xguide="Period length (hours)",guidefontsize=46)
plot!(xlimit=(0.0, 20.0))

plot!([4, 4], [0, 0.25], color=:grey)
plot!([8, 8], [0, 0.25], color=:grey)
plot!([16, 16], [0, 0.25], color=:grey)

# Plot example periodograms.
begin
    traj = rand(osc_trajs)
    x = 1/6:1/6:length(traj)/6
    p1 = plot(x, traj, xguide="Time (hours)", gridwidth=3,grid=true,xticks=[4:4:90...])
end

p2 = analyse_wpd_fit(traj)
plot(p1,p2,size=(4000,1600))

# Plot classification examples.
get_info = (t,s)->""
plot_classification_examples(exp_trajss[1]; info_f=get_info)
plot_classification_examples(exp_trajss[2]; info_f=get_info)
plot_classification_examples(exp_trajss[3]; info_f=get_info)
plot_classification_examples(exp_trajss[4]; info_f=get_info)
plot_classification_examples(exp_trajss[5]; info_f=get_info)

sim_traj_examples = rand(sim_trajss)
plot_classification_examples(sim_traj_examples; info_f=get_info)


### Extra ###
# Check distribution of base cases.
include("../Codebase/basic_analysis.jl")
include("../Codebase/file_management.jl")
include("../Codebase/utility.jl")
include("../Codebase/basic_plotting.jl")
include("../Codebase/models.jl")
include("../Codebase/model_simulation.jl")

@time sol = cle_sim([10.0,0.1, 2.0, 0.05, 3.0, 0.05], (0.0, 1000000.0))
traj = first.(sol.u)
# plot(traj)
wpd = welch_pgram(traj ./ mean(traj) .-1, 500)
plot(wpd.freq[2:end], wpd.power[2:end], label="Periodogram")

@. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
fit = curve_fit(model, wpd.freq[2:end], wpd.power[2:end], [0.5, 0.5])
plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.1); la=0.8,lw=7, label="Curve fitted to periodogram")

@time sol = cle_sim([10.0,10.0, 2.0, 0.05, 3.0, 0.01], (0.0, 1000.0))
traj = first.(sol.u)
# plot(traj)
wpd = welch_pgram(traj ./ mean(traj) .-1, 500)
plot(wpd.freq[2:end], wpd.power[2:end], label="Periodogram")

@. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
fit = curve_fit(model, wpd.freq[2:end], wpd.power[2:end], [0.5, 0.5])
plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.1); la=0.8,lw=7, label="Curve fitted to periodogram")


# Composite (full) distribution, (f, θ, σ, A, f0, s) = 4θ*(σ^2)/(4π*f^2 + θ^2) + A * exp(-(f-f0)^2/2(s^2))


### Old Stuff ###
osc_p_peaks = vcat(get_peak_periods.(osc_trajs)...)
hom_p_peaks = vcat(get_peak_periods.(hom_trajs)...)

density(osc_p_peaks; xguide="Period length (hours)",guidefontsize=46)
density!(hom_p_peaks)

osc_trajs = get_trajectory_class(exp_trajss[1], :oscillation)
hom_trajs = get_trajectory_class(exp_trajss[1], :homogeneous_activation)
osc_trajs = get_trajectory_class(exp_trajss[1], :oscillation)
periods = find_period.(osc_trajs)
osc_p_peaks = vcat(get_peak_periods.(osc_trajs)...)
density(osc_p_peaks; xguide="Period length (hours)",guidefontsize=46)


osc_trajs = get_trajectory_class(exp_trajss[5], :oscillation)
periods = find_period.(osc_trajs)
osc_p_peaks = vcat(get_peak_periods.(osc_trajs)...)
density!(osc_p_peaks; xguide="Period length (hours)",guidefontsize=46)
plot!(xlimit=(0.0, 20.0))

begin
    value_counts = countmap(periods)
    bar(collect(keys(value_counts)), collect(values(value_counts)))
end

plot_wpd(rand(osc_trajs))
plot_wpd(rand(hom_trajs))

traj = rand(osc_trajs)
analyse_wpd_fit(traj)

traj = rand(osc_trajs)
wpd = welch_pgram(traj ./ mean(traj) .-1, 200)
plot(wpd.freq,wpd.power; la=0.8,lw=7, label="Periodogram")
@. model(x, p) = p[1]*exp(-x*p[2])
fit = curve_fit(model, wpd.freq, wpd.power, [0.5, 0.5])
plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.1); la=0.8,lw=7, label="Exponential fitted to periodogram")

@. model(x, p) = 4*p[1]*(p[2]^2)/(4 * (pi^2) * (x^2) + (p[1]^2))
fit = curve_fit(model, wpd.freq[2:end], wpd.power[2:end], [1.0, 1.0])
plot!(x -> model(x, fit.param), 0.0, 0.5, xlimit=(0.0,0.1); la=0.8,lw=7, label="Exponential fitted to periodogram")

wpd.freq
minimum(wpd.power)

traj = rand(hom_trajs)
analyse_wpd_fit(traj)


osc_trajs = get_trajectory_class(exp_trajss[5], :oscillation)
osc_p_peaks = vcat(get_peak_periods.(osc_trajs)...)
density!(osc_p_peaks; xguide="Period length (hours)",guidefontsize=46)
plot!(xlimit=(0.0,25.0))

trajectory = exp_trajss[3][5][1]
plot(trajectory, xlimit=(1.0,100.0), xticks=1:6:100, grid=true, gridalpha=1)

trajectory = exp_trajss[3][5][1]
wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 150)
println(find_freq(wpd))
println("Period: $(round(1/(12*find_freq(wpd)), sigdigits=3))h")
plot(wpd.freq, wpd.power)

trajectory = exp_trajss[3][5][2]
wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 150)
println(find_freq(wpd))
println("Period: $(round(1/(12*find_freq(wpd)), sigdigits=3))h")
plot(wpd.freq, wpd.power)

trajectory = exp_trajss[3][5][3]
wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 250)
println(find_freq(wpd))
println("Period: $(round(1/(12*find_freq(wpd)), sigdigits=3))h")
plot(wpd.freq, wpd.power)

trajectory = exp_trajss[3][5][5]
wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 200)
println(find_freq(wpd))
println("Period: $(round(1/(12*find_freq(wpd)), sigdigits=3))h")
plot(wpd.freq, wpd.power)

trajectory = exp_trajss[3][5][6]
wpd = welch_pgram(trajectory ./ mean(trajectory) .-1, 200)
println(find_freq(wpd))
println("Period: $(round(1/(12*find_freq(wpd)), sigdigits=3))h")
plot(wpd.freq, wpd.power)

function find_freq(wpd)
    idx1 = findfirst(wpd.power[2:end] .> wpd.power[1:end-1])
    isnothing(idx1) && (return wpd.freq[1])
    idx2 = idx1 + argmax(wpd.power[idx1:end]) - 1
    isnothing(idx2) && (return wpd.freq[1])
    wpd.freq[idx2]
end

wpd.freq

### Old Code ###
using FFTW
using Plots
signal = exp_trajss[3][5][6] .- mean(exp_trajss[3][5][6])
fft_result = fft(signal)
power_spectrum = abs.(fft_result).^2
dominant_frequency_index = argmax(power_spectrum)
sample_rate = 1.0  # Assuming a unit sample rate (adjust if necessary)
dominant_frequency = dominant_frequency_index / length(signal) * sample_rate
frequency_axis = fftfreq(length(signal), 1 / sample_rate)
plot(frequency_axis, power_spectrum, xlabel="Frequency (Hz)", ylabel="Power Spectrum")
plot!(xlimit=(-0.02,0.02))

signal = exp_trajss[3][5][6]

fft_result = fft(signal)
power_spectrum = abs.(fft_result).^2
dominant_frequency_index = argmax(power_spectrum)
sample_rate = 1.0  # Assuming a unit sample rate (adjust if necessary)
dominant_frequency = dominant_frequency_index / length(signal) * sample_rate
frequency_axis = fftfreq(length(signal), 1 / sample_rate)
plot(frequency_axis, power_spectrum, xlabel="Frequency (Hz)", ylabel="Power Spectrum")


x = 0.0:0.05:314
signal = sin.(x)
wpd = welch_pgram(signal ./ mean(signal) .-1, 500)
println(wpd.freq[argmax(wpd.power)])
plot(wpd.freq, wpd.power)

idx = findfirst(wpd.power[2:end] .> wpd.power[1:end-1])
wpd.freq[idx+1]

2*0.008*pi/0.05

16*0.032

plot(x, signal)

using FFTW

# Assuming you have your sinusoidal signal in the 'signal' variable
signal = sin.(0.0:0.1:10.0) 
fft_result = fft(signal)
sample_rate = 1.0  # Sampling rate (adjust as needed)
frequency_bins = fftfreq(length(signal), 1 / sample_rate)
peak_index = argmax(abs.(fft_result))
peak_frequency = frequency_bins[peak_index]
plot(frequency_bins, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Amplitude Spectrum")
plot!()

using LombScargle

### Analysis ###
exp_trajss_all = load_all_trajectories()
exp_trajss = [[filter(e -> length(e) == 572, etjs) for etjs in exp_trajs] for exp_trajs in exp_trajss_all]
sim_trajss = load_all_simulations()


vals = exp_trajss[3][5][6]

vals = exp_trajss[3][5][6]
p1 = plot(vals)
plan = LombScargle.plan(1:length(vals), vals);
pgram = lombscargle(plan)
p2 = plot(freqpower(pgram)...)
plot(p1, p2, size=(3000,1200))


wpd = welch_pgram(vals)
plot(wpd.freq, wpd.power)

plan = LombScargle.plan(1:length(vals), vals);
pgram = lombscargle(plan)
plot(freqpower(pgram)...)



info_f = t -> ""
info_f = rescaled_wp_max
info_f = t -> percentile(t, 25)
plot_classification_examples(exp_trajss[1]; info_f=get_info)
plot_classification_examples(exp_trajss[2]; info_f=get_info)
plot_classification_examples(exp_trajss[3]; info_f=get_info)
plot_classification_examples(exp_trajss[4]; info_f=get_info)
plot_classification_examples(exp_trajss[5]; info_f=get_info)

plot_classification_examples(sim_trajs_1; info_f=get_info)
plot_classification_examples(sim_trajs_2; info_f=get_info)
plot_classification_examples(sim_trajs_3; info_f=get_info)
plot_classification_examples(sim_trajs_4; info_f=get_info)
plot_classification_examples(sim_trajs_5; info_f=get_info)

const sigv_data_dir = "Data/Experimental_data/Sigv_data"
DIR = "../Data/Experimental_data/SigV_data/Repeat_1/Corrected_one_step_1ugml_lysozyme_repeat_1.mat"
matform = read(matopen(DIR), "MY")
trajectories = [matform[:,i] for i in 1:size(matform,2)]
trajectories = [t[1:(any(isnan.(t)) ? findfirst(isnan.(t))-1 : length(t))] for t in trajectories]

plot(rand(trajectories[1:50])[50:end])

trajectories[50]

info_f=(t,s)->""
traj_set = exp_trajss[1]
set_data = get_set_data(traj_set)
class_transition_plot = plot_classification(traj_set; legend=:none)
classes = classify_trajectories(traj_set)
traj_set_classified = [traj_set[i][j] => classes[i][j] for i in 1:length(traj_set) for j in 1:length(traj_set[i])]
trajs_no_activation = sample_behaviour_type(traj_set_classified, :no_activation, 4)
trajs_stochastic_pulsing = sample_behaviour_type(traj_set_classified, :stochastic_pulsing, 4)
trajs_oscillation = sample_behaviour_type(traj_set_classified, :oscillation, 4)
trajs_homogeneous_activation = sample_behaviour_type(traj_set_classified, :homogeneous_activation, 4)
trajs_unknow = sample_behaviour_type(traj_set_classified, :unknown, 4)
trajs_additional = sample(traj_set_classified, 32-sum(length.([trajs_no_activation,trajs_stochastic_pulsing,trajs_oscillation,trajs_homogeneous_activation,trajs_unknow])); replace=false)
trajs_classified = vcat(trajs_no_activation,trajs_stochastic_pulsing,trajs_oscillation,trajs_homogeneous_activation,trajs_unknow, trajs_additional)
traj_class_plots = [plot_classified_trajectory(tc, set_data; info_f=info_f) for tc in trajs_classified]
plot(class_transition_plot, traj_class_plots..., layout=@layout[a{0.333w} grid(2,4); grid(4,6){0.65h}], size=(3000,2000))


trajs_oscillation = sample_behaviour_type(traj_set_classified, :oscillation, 100)
trajs = first.(trajs_oscillation)

trajectory = rand(trajs)
trajectory = trajectory .- mean(trajectory)
trajectory = trajectory ./ std(trajectory)
#trajectory = gaussian_smoothing(trajectory, 2.0)
wpgram = welch_pgram(trajectory, 100, 50)
p1 = plot(wpgram.freq, wpgram.power; yaxis=:log10,xlimit=(0.0,Inf))
p2 = plot(trajectory)
println(wpgram.freq[1])
plot(p1, p2, size=(3600,1200))

density(max_freq.(trajs))
histogram(max_freq.(trajs))
unique(max_freq.(trajs))

using ImageFiltering
function gaussian_smoothing(data::Vector{T}, σ::Float64) where T
    kernel = Kernel.gaussian(σ)
    return imfilter(data, kernel, "symmetric")
end


using StatsPlots
function max_freq(trajectory)
    trajectory = trajectory .- mean(trajectory)
    trajectory = trajectory ./ std(trajectory)
    wpgram = welch_pgram(trajectory .- mean(trajectory))
    arm = argmax(wpgram.power)
    wpgram.freq[arm]
end

# For a set of trajectories, plots their transition diagram, and 26 example trajectories and their classification.
function plot_classification_examples(traj_set::Vector{Vector{Vector{Float64}}}; info_f=(t,s)->"")
    set_data = get_set_data(traj_set)
    class_transition_plot = plot_classification(traj_set; legend=:none)
    classes = classify_trajectories(traj_set)
    traj_set_classified = [traj_set[i][j] => classes[i][j] for i in 1:length(traj_set) for j in 1:length(traj_set[i])]
    trajs_no_activation = sample_behaviour_type(traj_set_classified, :no_activation, 4)
    trajs_stochastic_pulsing = sample_behaviour_type(traj_set_classified, :stochastic_pulsing, 4)
    trajs_oscillation = sample_behaviour_type(traj_set_classified, :oscillation, 4)
    trajs_homogeneous_activation = sample_behaviour_type(traj_set_classified, :homogeneous_activation, 4)
    trajs_unknow = sample_behaviour_type(traj_set_classified, :unknown, 4)
    trajs_additional = sample(traj_set_classified, 32-sum(length.([trajs_no_activation,trajs_stochastic_pulsing,trajs_oscillation,trajs_homogeneous_activation,trajs_unknow])); replace=false)
    trajs_classified = vcat(trajs_no_activation,trajs_stochastic_pulsing,trajs_oscillation,trajs_homogeneous_activation,trajs_unknow, trajs_additional)
    traj_class_plots = [plot_classified_trajectory(tc, set_data; info_f=info_f) for tc in trajs_classified]
    plot(class_transition_plot, traj_class_plots..., layout=@layout[a{0.333w} grid(2,4); grid(4,6){0.65h}], size=(3000,2000),xticks=[],yticks=[])
end

plot_classification_examples(rand(sim_trajss); info_f=get_info)

sim_trajs_1 = rand(sim_trajss)
sim_trajs_2 = rand(sim_trajss)
sim_trajs_3 = rand(sim_trajss)
sim_trajs_4 = rand(sim_trajss)
sim_trajs_5 = rand(sim_trajss)

set_data = get_set_data(sim_trajs_1)
(set_data[:min_val] + set_data[:max_val])/5

get_info(trajectory, set_data) = ""
function get_info(trajectory, set_data)    
    "$(round(percentile(trajectory, 80),sigdigits=3))  <  $(round(set_data[:max_val]/5,sigdigits=3))
     $(round(maximum(trajectory),sigdigits=3))  <  $(round(set_data[:max_val]/5,sigdigits=3))"
end

savefig("sim2.png")


plot( 1:5 )

img = rand(10, 10)
p1 = heatmap(img, color=:grays, legend=false)
plot!( -5:8, (-5:8).^2, inset = (1, bbox(0.1,0.0,0.4,0.4)), subplot = 2)


using Plots

x = 1:10
y = sin.(x)

# Main heatmap plot
img = rand(10, 10)
p1 = heatmap(img, color=:grays, legend=false)

plot!(x, y, linewidth=2, color=:red, legend=false,
background_color=:transparent, 
subplot_background_color=:transparent,
grid=false, inset = (1, bbox(0.1,0.0,0.4,0.4)), subplot = 2)

# Overlay plot with transparent background
p2 = plot(x, y, linewidth=2, color=:red, legend=false,
          background_color=:transparent, 
          subplot_background_color=:transparent,
          grid=false)

# Overlay p2 over p1
plot!(p1, inset = (p2, bbox(0.0, 0.0, 1.0, 1.0)), subplot = 2)


### Load matlab SigV trajectories ###

DIR = "../Data/Experimental_data/SigV_data/Repeat_1/Corrected_one_step_1ugml_lysozyme_repeat_1.mat"
data = read(matopen(DIR), "MY")
trajectories_unfiltered = [data[:, i] for i in 1:size(data, 2)]

function filter_Nans(trajectory)
    idx = findfirst(isnan, trajectory)
    if idx === nothing
        return trajectory
    else
        return trajectory[1:idx-1]
    end
end
trajectories = filter(traj -> !isempty(traj), filter_Nans.(trajectories_unfiltered))

dir = "../Data/Experimental_data/SigV_data/Repeat_1/"
for file in readdir(dir)
    file_path = dir*file
    data = read(matopen(file_path), "MY")
    trajectories_unfiltered = [data[:, i] for i in 1:size(data, 2)]
    
    function filter_Nans(trajectory)
        idx = findfirst(isnan, trajectory)
        if idx === nothing
            return trajectory
        else
            return trajectory[1:idx-1]
        end
    end
    trajectories = filter(traj -> !isempty(traj), filter_Nans.(trajectories_unfiltered))
    serialize(replace(file_path,".mat" => ".jls"),trajectories)
end

readdir("../Data/Experimental_data/SigV_data/Repeat_1/")

data = deserialize("../Data/Experimental_data/SigV_data/Repeat_2/Corrected_one_step_1ugml_lysozyme_repeat_2.jls")
plot(data, color=1, la=0.6)

data = deserialize("../Data/Experimental_data/SigV_data/Repeat_2/Corrected_one_step_2ugml_lysozyme_repeat_2.jls")
plot!(data, color=2, la=0.6)
data = deserialize("../Data/Experimental_data/SigV_data/Repeat_2/Corrected_one_step_4ugml_lysozyme_repeat_2.jls")
plot!(data, color=3, la=0.6)
plot!(xlimit=(20,30),xticks=1:40)

mf = matopen(DIR)
names(mf) |> println
println(read(mf, "time_zero")[1:200,1])
read(mf, "MY")