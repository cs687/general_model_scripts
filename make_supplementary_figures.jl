####################
### Preparations ###
####################

### Generates Occurrence Counts ####

### Fetch Functions ###
include("../Codebase/basic_analysis.jl")
include("../Codebase/file_management.jl")
include("../Codebase/utility.jl")
include("../Codebase/basic_plotting.jl")
include("../Codebase/models.jl")
include("../Codebase/model_simulation.jl")
include("../Codebase/behaviour_grid_loading.jl")
include("../Codebase/behaviour_grid_plotting.jl")
include("../Codebase/make_occurence_analysis.jl")

### Behaviour Map Parameters ###
S_grid = 10 .^(range(-1,stop=2,length=300))
D_grid = 10 .^(range(-1,stop=2,length=300))
τ_grid = [0.1,0.15,0.20,0.30,0.50,0.75,1.0,1.5,2.0,3.0,5.0,7.50,10.0,15.0,20.0,30.0,50.0,75.0,100.0]
v0_grid = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.20]
n_grid = [2.0,3.0,4.0]
η_grid = [0.001,0.002,0.005,0.01,0.02,0.05,0.1]

target_folder = "BehaviourMap"
dataset = DataSet(target_folder,S_grid,D_grid,τ_grid,v0_grid,n_grid,η_grid);
behaviour_occurences = count_occurences(dataset);


### Plot Settings ###

# Default values.
pythonplot();
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=16, yguidefontsize=16, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);

###############
### Figures ###
###############

### Supplementary Figure 1, Nullclines ###
let 
    # Subfigure 1.
    parameters_1 = [1.65, 0.2, 10.0, 0.1, 2, 0.0]
    plot_nullcline_set(parameters_1, ymax=1, lw=8)
    qLx = 0.06; qLy = 0.04;; qW = 2;
    (x,y)=(0.45,0.6); (d1,d2)=(1,-1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.3,0.85); (d1,d2)=(-1,-1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.1,0.6); (d1,d2)=(1,-1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.8,0.4); (d1,d2)=(-1,1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.45,0.25); (d1,d2)=(1,1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.25,0.1); (d1,d2)=(-1,1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    (x,y)=(0.1,0.0225); (d1,d2)=(1,1); scatter!((x,y),color=:black,label=""); quiver!([x,x],[y,y],quiver=([0.0,d1*qLy],[d2*qLx,0.0]),color=:black,size=(600,400),lw=qW)
    p1 = plot!(xticks=[],ytick=[],left_margin=4mm)

    # Subfigure 2.
    parameters_2 = [1.65, 0.2, 10.0, 0.1, 2, 0.0]
    p2 = plot_nullcline_sets(parameters_2, :S, 1.6:0.05:1.8; xmax=1.0, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)

    # Subfigure 3.
    parameters_3 = [1.7, 0.2, 10.0, 0.1, 2, 0.0]
    p3 = plot_nullcline_sets(parameters_3, :D, 0.1:0.05:0.3, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)

    # Subfigures 4, 5.
    parameters_4 = [1.65, 0.2, 10.0, 0.1, 2, 0.0]
    p4 = plot_nullcline_sets(parameters_4, :S, 1.6:0.05:1.8, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)
    parameters_5 = [1.65, 0.2, 10.0, 0.1, 3, 0.0]
    p5 = plot_nullcline_sets(parameters_5, :S, 1.6:0.05:1.8, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)

    # Subfigures 6, 7, 8.
    parameters_6 = [1.0, 0.2, 10.0, 0.1, 2, 0.0]
    p6 = plot_nullcline_sets(parameters_6, :v0, 0.05:0.05:0.25; xmax=1.5, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)
    parameters_7 = [2.0, 0.2, 10.0, 0.1, 2, 0.0]
    p7 = plot_nullcline_sets(parameters_7, :v0, 0.05:0.05:0.25; xmax=1.5, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)
    parameters_8 = [3.0, 0.2, 10.0, 0.1, 2, 0.0]
    p8 = plot_nullcline_sets(parameters_8, :v0, 0.05:0.05:0.25; xmax=1.5, ymax=5., legendfontsize=12, xticks=[], yticks=[], left_margin=4mm, lw=5)
    asfasf = 5
    nullcline_figure = plot(p1,p2,p3,p6,p7,p8,p4,p5,size=(1600,800),layout=@layout [a b c; d{0.7h} e f g h])
    make_figure(nullcline_figure, "Nullclines"; supp=true)
end

#end

### Supplementary Figure 2, Bifurcation Diagrams ###
let 
    vert_line_h = 1.25

    parameters1 = [0.0, 0.25, 10.0, 0.1, 3, 0.]
    S_span1 = (1.0, 4.0)
    S_vals1 = [1.25, 1.75, 2.25, 3.0]
    bif1 = plot_bifurcation_diagram(parameters1, S_span1, :S; make_labels=true, lw=9, la=0.8)
    foreach(i -> plot_vertical_line!(S_vals1[i], vert_line_h, color=red_scale[i], lw=5, linestyle=:dash, la=0.7), 1:length(S_vals1))
    p1 = plot!(xlimit=S_span1,ylimit=(0.,vert_line_h),xticks=[],yticks=[])
    p2 = plot_nullcline_sets(parameters1, :S, S_vals1; make_labels=true, colors=red_scale,legendfontsize=12,xticks=[],yticks=[],xmax=1.25,ymax=10.0,lw=7,la=0.8)
    
    
    parameters2 = [0.,5.,10.0,0.1,3,0.]
    S_span2 = (1.,12.)
    S_vals2 = [2,3.5,5,8]
    bif2 = plot_bifurcation_diagram(parameters2, S_span2, :S; make_labels=true, lw=9, la=0.8)
    foreach(i -> plot_vertical_line!(S_vals2[i],vert_line_h,color=red_scale[i],lw=5, linestyle=:dash, la=0.7), 1:length(S_vals2))
    p3 = plot!(xlimit=S_span2,ylimit=(0.,vert_line_h),xticks=[],yticks=[])
    p4 = plot_nullcline_sets(parameters2, :S, S_vals2; make_labels=true, colors=red_scale,legendfontsize=12,xticks=[],yticks=[],xmax=1.25,ymax=1.25,lw=7,la=0.8)
    
    S_vals3 = [1.75,2.5,3,3.85]
    S_span3 = (1.,4.)
    parameters3 = [0.,2.8,20.0,0.1,3,0.]
    bif3 = plot_bifurcation_diagram(parameters3, S_span3, :S; make_labels=true, lw=9, la=0.8)
    foreach(i -> plot_vertical_line!(S_vals3[i],vert_line_h,color=red_scale[i],lw=5, linestyle=:dash, la=0.7), 1:length(S_vals3))
    p5 = plot!(xlimit=S_span3,ylimit=(0.,vert_line_h),xticks=[],yticks=[])
    p6 = plot_nullcline_sets(parameters3, :S, S_vals3; make_labels=true, colors=red_scale,legendfontsize=12,xticks=[],yticks=[],xmax=1.25,ymax=1.25,lw=7,la=0.8);
    
    p7 = plot_bifurcation_diagrams([1.,.25,100.0,0.1,3,0],(1.,12.),:S,[0.25,1.,2.,4.],:D; legendfontsize=12, lw=6, colours_s = bluegreen_scale,colours_u=red_scale)
    
    p8 = plot_bifurcation_diagrams([0.,0.25,10.0,0.01,3,0.],(1.,8.),:S,[0.025,0.05,0.1,0.2],:v0; legendfontsize=12, lw=6, colours_s = bluegreen_scale,colours_u=red_scale)
    
    p9 = plot_bifurcation_diagrams([0.,0.25,10.0,0.01,3,0.],(1.,12.),:S,[1.0,2.0,3.0],:n; legendfontsize=12, lw=6, colours_s = bluegreen_scale,colours_u=red_scale)
    
    p101 = plot_bifurcation_diagram([1., 4.0, 0.1, 0.1, 3, 0.0], (1.,12.), :S, title="tau=0.1", lw=6)
    p102 = plot_bifurcation_diagram([1., 4.0, 0.5, 0.1, 3, 0.0], (1.,12.), :S, title="tau=0.5", lw=6)
    p103 = plot_bifurcation_diagram([1., 4.0, 2.0, 0.1, 3, 0.0], (1.,12.), :S, title="tau=2.0", lw=6)
    p104 = plot_bifurcation_diagram([1., 4.0, 10.0, 0.1, 3, 0.0], (1.,12.), :S, title="tau=10.0", lw=6)
    p10 = plot(p101,p102,p103,p104,size=(400,500),layout=(4,1),title="",xguide="",yguide="",ylimit=(0.,1.25))
    
    bifurcation_diagrams_figure = plot(p1,p5,p3,p2,p6,p4,p7,p8,p9,p10,size=(2000,1000),xticks=[],yticks=[],bottom_margin=5mm,left_margin=5mm,layout=@layout [grid(1,3); grid(1,3){0.28h}; grid(1,4)])
    make_figure(bifurcation_diagrams_figure, "Bifurcation_diagrams"; supp=true)
end


### Supplementary Figure 3, Behaviour Map with Bifurcation Diagrams ###
let
    # Upper half (uses `gr`).
    gr()
    params = [2.0,0.1,3,0.1]
    bg = BehaviourGrid(params, dataset)
    plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="D", yguide="S")
    plot!([40,40],[0.5,200.],color=:black,l2=1,label="", lw=3); plot!([40,40],[0.5,200.],color=:white,l2=1,label="", lw=1)
    plot!([100,100],[0.5,200.],color=:black,l2=1,label="", lw=3); plot!([100,100],[0.5,200.],color=:white,l2=1,label="", lw=1)
    plot!([170,170],[0.5,200.],color=:black,l2=1,label="", lw=3); plot!([170,170],[0.5,200.],color=:white,l2=1,label="", lw=1)
    plot!([240,240],[0.5,200.],color=:black,l2=1,label="", lw=3); plot!([240,240],[0.5,200.],color=:white,l2=1,label="", lw=1)
    p1 = plot!(xlimit=(-Inf,Inf),ylimit=(-Inf,Inf),left_margin=0mm,bottom_margin=0mm)

    behaviour_grid_and_bifurcation_figure_top = plot(p0_gr,p1,p0_gr,size=(1500,600),layout=@layout [a{0.25w} b c{0.25w}])
    make_figure(behaviour_grid_and_bifurcation_figure_top, "Behaviour_grid_and_bfurcation_diagrams"; supp=true, tag="top")

    # Upper half (uses `pythonplot`).
    pythonplot()
    p2 = plot_bifurcation_diagram([0.,0.2,params...],(dataset.S_grid[100],100.),:S,title="D=0.2", xscale=:log10,lw=7,lsU=:dot,la=0.8,dsmax=0.001,maxSteps=10000000)
    p3 = plot_bifurcation_diagram([0.,1.0,params...],(dataset.S_grid[100],100.),:S,title="D=1.0", xscale=:log10,lw=7,lsU=:dot,la=0.8,dsmax=0.001,maxSteps=10000000)
    p4 = plot_bifurcation_diagram([0.,5.,params...],(dataset.S_grid[100],100.),:S,title="D=5.0", xscale=:log10,lw=7,lsU=:dot,la=0.8)
    p5 = plot_bifurcation_diagram([0.,25.,params...],(dataset.S_grid[100],100.),:S,title="D=25.0", xscale=:log10,lw=7,lsU=:dot,la=0.8)

    plot(p2,p3,p4,p5,layout=(1,4),size=(1600,350),bottom_margin=0mm,left_margin=0mm)

    behaviour_grid_and_bifurcation_figure_bottom = plot(p2,p3,p4,p5,size=(1500,300),layout=(1,4))
    make_figure(behaviour_grid_and_bifurcation_figure_bottom, "Behaviour_grid_and_bfurcation_diagrams"; supp=true, tag="bottom")
end

### Supplementary Figure 4, Behaviour Map Examples ###
let
    gr()
    p1 = [2.0,0.1,3,0.1]; τs = [τ_grid[1],τ_grid[5],τ_grid[6],τ_grid[9],τ_grid[15]]
    bgs = map(τ -> BehaviourGrid(setindex!(copy(p1),τ,1),dataset), τs)
    bg_plots = map(bg -> plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="",yguide="",title="τ=$(bg.params[1])"), bgs)
    p1 = plot(bg_plots...,size=(1600,300),layout=(1,5))

    p2 = [2.0,0.1,3,0.1]; v₀s = [v0_grid[2],v0_grid[4],v0_grid[6],v0_grid[7],v0_grid[8]]
    bgs = map(v₀ -> BehaviourGrid(setindex!(copy(p2),v₀,2),dataset), v₀s)
    bg_plots = map(bg -> plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="",yguide="",title="v₀=$(bg.params[2])"), bgs)
    p2 = plot(bg_plots...,size=(1400,300),layout=(1,5))

    p3 = [2.0,0.1,3,0.1]; ηs = [η_grid[1],η_grid[3],η_grid[5],η_grid[6],η_grid[7]]
    bgs = map(η -> BehaviourGrid(setindex!(copy(p3),η,4),dataset), ηs)
    bg_plots = map(bg -> plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="",yguide="",title="η=$(bg.params[4])"), bgs)
    p3 = plot(bg_plots...,size=(1400,300),layout=(1,5))

    p4 = [2.0,0.1,3,0.1]; ns = [n_grid[1],n_grid[2],n_grid[3]]
    bgs = map(n -> BehaviourGrid(setindex!(copy(p4),n,3),dataset), ns)
    bg_plots = map(bg -> plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="",yguide="",title="n=$(bg.params[3])"), bgs)
    p4 = plot(p0_gr,bg_plots...,p0_gr,size=(1400,300),layout=(1,5))

    behaviour_map_examples_figure = plot(p1,p2,p3,p4,size=(2000,1200),xticks=[],yticks=[],layout=@layout [a ; b; c ; d])
    make_figure(behaviour_map_examples_figure, "Behaviour_map_examples"; supp=true)
end


### Supplementary Figure 5, τ-η Maps ###
let 
    gr()
    params = [2.0,0.1,3,0.05];
    bgs = BehaviourGrids(params,dataset;τ_values=[τ_grid[1], τ_grid[3], τ_grid[4], τ_grid[5], τ_grid[7], τ_grid[11], τ_grid[15]],η_values=[η_grid[1],η_grid[3],η_grid[4],η_grid[5],η_grid[7]]);
    τ_η_maps_figure = plot_behaviour_grids(bgs;plot_size=(1860,1685),start_s_slice=100,idx_axis=false)    
    τ_η_maps_figure = plot!(τ_η_maps_figure; left_margin=3mm, bottom_margin=3mm)
    make_figure(τ_η_maps_figure, "Tau_eta_maps"; supp=true)
end


### Supplementary Figure 6, τ-v0 Maps ###
let 
    gr()
    params = [2.0,0.1,3,0.05];
    bgs = BehaviourGrids(params,dataset;τ_values=[τ_grid[1], τ_grid[3], τ_grid[4], τ_grid[5], τ_grid[7], τ_grid[11], τ_grid[15]],v0_values=[v0_grid[1],v0_grid[3],v0_grid[4],v0_grid[6],v0_grid[8]]);
    τ_v0_maps_figure = plot_behaviour_grids(bgs;plot_size=(1860,1685),start_s_slice=100,idx_axis=false)
    τ_v0_maps_figure = plot!(τ_v0_maps_figure; left_margin=3mm, bottom_margin=3mm)
    make_figure(τ_v0_maps_figure, "Tau_v0_maps"; supp=true)
end


### Supplementary Figure 7, v0-η Maps ###
let 
    gr()
    params = [2.0,0.1,3,0.05];
    bgs = BehaviourGrids(params,dataset;v0_values=[v0_grid[1],v0_grid[2],v0_grid[4],v0_grid[5],v0_grid[6],v0_grid[7],v0_grid[8]],η_values=[η_grid[1],η_grid[3],η_grid[4],η_grid[5],η_grid[7]]);
    v0_η_maps_figure = plot_behaviour_grids(bgs;plot_size=(1860,1685),start_s_slice=100,idx_axis=false,transpose=true)
    v0_η_maps_figure = plot!(v0_η_maps_figure; left_margin=3mm, bottom_margin=3mm)
    make_figure(v0_η_maps_figure, "v0_eta_maps"; supp=true)
end


### Supplementary Figure 8, τ-n Maps ###
let 
    gr()
    params = [2.0,0.1,3,0.05];
    bgs = BehaviourGrids(params,dataset;τ_values=[τ_grid[1], τ_grid[3], τ_grid[4], τ_grid[5], τ_grid[7], τ_grid[11], τ_grid[15]],n_values=[n_grid[1],n_grid[2],n_grid[3]]);
    τ_n_maps_figure = plot_behaviour_grids(bgs;plot_size=(1115,1680),start_s_slice=100,idx_axis=false)
    τ_n_maps_figure = plot!(τ_n_maps_figure; left_margin=3mm, bottom_margin=3mm)
    make_figure(τ_n_maps_figure, "Tau_n_maps"; supp=true)
end


### Supplementary Figures 9-19, Behaviour descriptions ###
begin 
    σmax = 1.25; Amax = 1.25;
    time_pts = (-100.0, 0.0, 500.0)
    n = 3;
    la = 0.45
    xguide=""; yguide="";
end

# No activation.
let 
    pythonplot()
    parameters_no_activation = [1.0,  10.0, 10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_no_activation, time_pts, n, la=la, ymax=σmax, lw=3, yguide="Concentration (σ)")
    p2,p3 = plot_activation_time_n_pspace(parameters_no_activation, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
    make_figure(plot_no_activation, "Behaviour_No_activation"; supp=true)
end

include("../Codebase/make_occurence_analysis.jl")
# Stochastic pulsing.
let 
    pythonplot()
    parameters_stochastic_pulsing = [2.8,  10.0, 10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_stochastic_pulsing, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_stochastic_pulsing, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:stochastic_pulsing, behaviour_occurences, dataset)    
    plot_stochastic_pulsing = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_stochastic_pulsing, "Behaviour_Stochastic_pulsing"; supp=true)
end

# Oscillation.
let 
    pythonplot()
    parameters_oscillation = [10.0, 10.0, 10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_oscillation, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_oscillation, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:oscillation, behaviour_occurences, dataset)    
    plot_oscillation = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_oscillation, "Behaviour_Oscillation"; supp=true)

end

# Stochastic anti-pulsing.
let 
    pythonplot()
    parameters_stochastic_anti_pulsing = [15.5, 10.0, 10.0, 0.05,  3, 0.025]
    p1 = plot_activations(parameters_stochastic_anti_pulsing, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_stochastic_anti_pulsing, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:stochastic_anti_pulsing, behaviour_occurences, dataset)    
    plot_stochastic_anti_pulsing = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_stochastic_anti_pulsing, "Behaviour_Stochastic_anti_pulsing"; supp=true)
end

# Heterogeneous activation.
let 
    pythonplot()
    parameters_heterogeneous_activation = [2.25, 0.1,  10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_heterogeneous_activation, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_heterogeneous_activation, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:heterogeneous_activation, behaviour_occurences, dataset)    
    plot_heterogeneous_activation = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_heterogeneous_activation, "Behaviour_Heterogeneous_activation"; supp=true)
end

# Homogeneous activation.
let 
    pythonplot()
    parameters_homogeneous_activation = [5.0,  0.1,  10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_homogeneous_activation, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_homogeneous_activation, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:homogeneous_activation, behaviour_occurences, dataset)    
    plot_homogeneous_activation = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_homogeneous_activation, "Behaviour_Homogeneous_activation"; supp=true)
end

# Stochastic switching.
let 
    pythonplot()
    parameters_stochastic_switching = [2.25, 1.2,  10.0, 0.05,  2, 0.025]
    p1 = plot_activations(parameters_stochastic_switching, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_stochastic_switching, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:stochastic_switching, behaviour_occurences, dataset)    
    plot_stochastic_switching = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_stochastic_switching, "Behaviour_Stochastic_switching"; supp=true)
end

# Stable bistability.
let 
    pythonplot()
    parameters_stable_bistability = [2.1,  0.1,  10.0, 0.05,  2, 0.025]
    plot_activations(parameters_stable_bistability, time_pts, 1; colors=blue_scale[1], la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    plot_simulation!(parameters_stable_bistability, (0.0, time_pts[3]), idxs=[1], color=blue_scale[3], ymax=Amax, la=0.5, u0=[0.8,0.8,0.8,0.8], lw=3, xticks=:auto, yticks=:auto)
    p1 = plot!(xlimit=(time_pts[1],time_pts[3]))

    sim1 = cle_sim(parameters_stable_bistability,(time_pts[1],time_pts[3])); 
    sim2 = cle_sim(parameters_stable_bistability,(0.0,time_pts[3]),u0=[0.8,0.8,0.8,0.8]);
    plot(sim1,vars=[1,4],color=[:deepskyblue :orange],la=[0.9 0.7],lw=[4 3],label=""); 
    plot!(sim2,vars=[1,4],ylimit=(0.,σmax),color=[:darkblue :darkorange3],la=[0.9 0.7],lw=[4 3],label="")
    p2 = plot_vertical_line!(0,σmax,lw=7,la=[0.8 0.75],left_margin=5mm,xlimit=(time_pts[1],time_pts[3]),xguide="Time",xticks=[],yticks=[])
    p2 = plot!(p2, yguide="", yticks=[])
    plot(sim1,vars=(1,4),la=0.8,lw=4,color=:deepskyblue,label=""); 
    plot!(sim2,vars=(1,4),xlimit=(0.,σmax),ylimit=(0.,Amax),la=0.8,lw=4,color=:darkblue,label="")
    p3 = plot_nullcline_sets!(parameters_stable_bistability,:S,[parameters_stable_bistability[4],parameters_stable_bistability[1]],color1=:green,colors2=[:pink, :red],lw=5,set_label=false,xmax=σmax,ymax=Amax,title="",left_margin=5mm,xticks=[],yticks=[])

    p4 = plot_occurences(:stable_bistability, behaviour_occurences, dataset)
    plot_stable_bistability = plot(p1,p2,p3,p4...,size=(1500,800),bottom_margin=5mm,layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_stable_bistability, "Behaviour_Stable_bistability"; supp=true)
end

# Single response-pulse.
let 
    pythonplot()
    parameters_single_response_pulse = [2.75, 10.0, 10.0, 0.05,  2, 0.01]
    p1 = plot_activations(parameters_single_response_pulse, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_single_response_pulse, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:single_response_pulse, behaviour_occurences, dataset)    
    plot_single_response_pulse = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_single_response_pulse, "Behaviour_Single_response_pulse"; supp=true)
end

# Intermediary activation.
let 
    pythonplot()
    parameters_intermediary_activation = [8.0,  10.0, 0.1,  0.05,  2, 0.025]
    p1 = plot_activations(parameters_intermediary_activation, time_pts, n, la=la, ymax=σmax, lw=3, xticks=:auto, yticks=:auto)
    p2,p3 = plot_activation_time_n_pspace(parameters_intermediary_activation, time_pts,ymax_t=σmax, xmax_ps=σmax, ymax_ps=Amax);
    p2 = plot!(p2, yguide="", yticks=[])
    plot_no_activation = plot(p1,p2,p3,size=(1500,480),bottom_margin=7mm,left_margin=9mm,layout=@layout [a b c])
   
    p4 = plot_occurences(:homogeneous_intermediate_activation, behaviour_occurences, dataset)    
    plot_intermediary_activation = plot(p1,p2,p3,p4...,size=(1500,800),layout=@layout [a b c; d{0.4h} e f g])
    make_figure(plot_intermediary_activation, "Behaviour_Intermediary_activation"; supp=true)
end


### Supplementary Figure 20, Gillespie simulations ###
@time let 
    time_pts = (-12000.0, 0.0, 50000.0)
    ymax = 250
    la=0.9

    # [v0, v, I, S, D, K, n, d, τ]
    parameters_no_activation =              [0.1,  5.0,  1.0, 0.5,   150.0, 100.0, 2, 0.01,  500.0]
    parameters_stochastic_pulsing =         [0.1,  2.0,  1.0, 15.0,  150.0, 10.0,  3, 0.01,  1500.0];
    parameters_oscillation =                [0.1,  2.5,  1.0, 120.0, 150.0, 10.0,  3, 0.01,  1000.0];
    parameters_stochastic_anti_pulsing =    [0.1,  2.0,  1.0, 240.0, 150.0, 10.0,  3, 0.01,  500.0];
    parameters_heterogeneous_activation =   [0.1,  2.0,  1.0, 0.25,  0.01,  15.0,  3, 0.01,  500.0];
    parameters_homogeneous_activation =     [0.1,  2.0,  1.0, 1.0,   0.01,  15.0,  3, 0.01,  500.0];
    parameters_stochastic_switching =       [0.1,  1.2,  1.0, 0.02,  0.01,  1.0,   3, 0.01,  100.0];
    parameters_stable_bistability =         [0.05, 2.0,  1.0, 0.08,  0.01,  5.0,   3, 0.01,  500.0];
    parameters_single_pulse_response =      [0.25, 5.0,  1.0, 38.0,  150.0, 10.0,  3, 0.001, 1000.0];
    parameters_intermediary_activation =    [0.1,  3.5,  1.0, 100.0, 150.0, 10.0,  3, 0.01,  5.0];

    no_activation_ssa_plot = plot_activation(parameters_no_activation, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    stochastic_pulsing_ssa_plot = plot_activation(parameters_stochastic_pulsing, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    oscillation_ssa_plot = plot_activation(parameters_oscillation, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    stochastic_anti_pulsing_ssa_plot = plot_activation(parameters_stochastic_anti_pulsing, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    heterogeneous_activation_ssa_plot = plot_activations(parameters_heterogeneous_activation, time_pts, 3; gillespie=true, la=la, ymax=ymax,idxs=1)
    homogeneous_activation_ssa_plot = plot_activation(parameters_homogeneous_activation, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    stochastic_switching_ssa_plot = plot_activation(parameters_stochastic_switching, 5.0 .* time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    plot_activation(parameters_stable_bistability, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)
    plot_simulation!(parameters_stable_bistability, (time_pts[2],time_pts[3]); gillespie=true, la=la, idxs=1, ymax=ymax,color=blue_scale[3], u0=[200,200,200,200])
    stable_bistability_ssa_plot = plot!(xlimit=(time_pts[1],time_pts[3])) 
    single_pulse_response_ssa_plot = plot_activation(parameters_single_pulse_response, 5.0 .* time_pts; gillespie=true, la=la, ymax=5000,idxs=1)
    intermediary_activation_ssa_plot = plot_activation(parameters_intermediary_activation, time_pts; gillespie=true, la=la, ymax=ymax,idxs=1)

    gillespie_behaviour_recreation_figure = plot(no_activation_ssa_plot,stochastic_pulsing_ssa_plot,oscillation_ssa_plot,stochastic_anti_pulsing_ssa_plot,heterogeneous_activation_ssa_plot,homogeneous_activation_ssa_plot,stochastic_switching_ssa_plot,stable_bistability_ssa_plot,single_pulse_response_ssa_plot,intermediary_activation_ssa_plot,yguide="Copy numbers",xlimit=(-10000.0, 50000.0),layout=(2,5),size=(3500,800))
    make_figure(gillespie_behaviour_recreation_figure, "Gillespie_behaviour_recreation"; supp=true)
end

# Grid classification for Gillespie interpretation.
let 
    gr()
    include("../Codebase/Functions_specific/gillespie_behaviour_grids.jl")
    S_grid = 10 .^(range(-2,stop=2,length=300))
    D_grid = 10 .^(range(-2,stop=2,length=300))

    gillespie_grid_classification = deserialize("../Data/Simulated_data/GillespieMap/gillespie_grid_classification.jls")
    # Ticks removed since these does not necessarily correpsond with what is correct (this case is a bit special).
    behaviour_grid_plot = plot!(plot_behaviour_grid(gillespie_grid_classification, S_grid, D_grid; idx_axis=false, xguide="", yguide="", xticks=[], yticks=[]),size=(600,400))   
    make_figure(behaviour_grid_plot, "Gillespie_classification_map"; supp=true)
end

### Supplementary Figure 11, Single A in Time Delay Chain ###
let 
    include("../Codebase/Functions_specific/intermediary_numbers_investigation.jl")
    la=0.9

    # Plots simulations with various number of intermediaries. 
    parameters = [5.,5.,20.0,0.1,4.,0.0];
    plot_oscillation_simulation(parameters, model1, 210., label="1 Intermediary")
    plot_oscillation_simulation!(parameters, model2, 210., color=2,label="2 Intermediaries")
    plot_oscillation_simulation!(parameters, model3, 210., color=3,label="3 Intermediaries")
    plot_oscillation_simulation!(parameters, model4, 210., color=4,label="4 Intermediaries");
    p1 = plot!(legend=:bottomright, yguide="σ", ylimit=(0.0,1.25))

    # Analyses period length as function of number of inetermediaries.
    τ_vals = 10 .^ range(log10(10.0),stop=log10(1000.0),length=20)
    int_nbrs = 1:4
    periods = map(i -> map(τ -> meassure_period(parameters,models[i],τ=τ), τ_vals), int_nbrs)
    plot(); foreach(idx -> plot!(τ_vals,periods[idx],label="$(idx) Intermediaries",la=la, lw=4), 1:length(periods))
    p2 = plot!(xlimit=(τ_vals[1],τ_vals[end]),legend=:bottomright, xaxis=:log10,xguide="τ",yguide="Period length", ylimit=(0.0,2500.0))

    intermediaries_number_analysis_figure = plot(p1,p2,size=(1500,350))
    make_figure(intermediaries_number_analysis_figure, "Intermediaries_number_analysis"; supp=true)
end


### Supplementary Figure 12, SDDE ###
let 
    include("../Codebase/Functions_specific/model_sdde_simulations.jl")

    t_pts = (-100.0, 0.0, 500.0)
    ymax = 1.25;
    
    parameters_sdde_no_activation =             [1.0,  10.0, 25.0, 0.05,  2, 0.025]
    parameters_sdde_stochastic_pulsing =        [2.5,  10.0, 25.0, 0.05,  2, 0.05]
    parameters_sdde_oscillation =               [5.0,  10.0, 25.0, 0.05,  2, 0.025]
    parameters_sdde_stochastic_anti_pulsing =   [20.0, 10.0, 25.0, 0.025,  2, 0.1]
    parameters_sdde_heterogeneous_activation =  [2.25, 0.1,  25.0, 0.05,  2, 0.025]
    parameters_sdde_homogeneous_activation =    [3.0,  0.1,  25.0, 0.05,  2, 0.025]
    parameters_sdde_stochastic_switching =      [2.22,  0.9, 25.0, 0.05,  2, 0.05]
    parameters_sdde_stable_bistability =        [2.0,  0.1,  25.0, 0.05,  2, 0.025]
    parameters_sdde_single_pulse_response =     [3.6,  10.0, 25.0, 0.25,  3, 0.02]
    parameters_sdde_intermediary_activation =   [7.5,  10.0, 1.0,  0.05,  2, 0.025]
    
    plot_sdde_no_activation = plot_sdde_activation(parameters_sdde_no_activation, t_pts, ymax=ymax)
    plot_sdde_stochastic_pulsing = plot_sdde_activation(parameters_sdde_stochastic_pulsing, t_pts, ymax=ymax)
    plot_sdde_oscillation = plot_sdde_activation(parameters_sdde_oscillation, t_pts, ymax=ymax)
    plot_sdde_stochastic_anti_pulsing = plot_sdde_activation(parameters_sdde_stochastic_anti_pulsing, t_pts, ymax=ymax)
    plot_sdde_activation(parameters_sdde_heterogeneous_activation, t_pts, ymax=ymax, color=blue_scale[1], la=0.7)
    plot_sdde_activation!(parameters_sdde_heterogeneous_activation, t_pts, ymax=ymax, color=blue_scale[2], la=0.6)
    plot_sdde_heterogeneous_activation = plot_sdde_activation!(parameters_sdde_heterogeneous_activation, t_pts, ymax=ymax, color=blue_scale[3], la=0.5)
    plot_sdde_homogeneous_activation = plot_sdde_activation(parameters_sdde_homogeneous_activation, t_pts, ymax=ymax)
    plot_sdde_stochastic_switching = plot_sdde_activation(parameters_sdde_stochastic_switching, t_pts, ymax=ymax)
    plot_sdde_activation(parameters_sdde_stable_bistability, (0.0, 0.0001, 500.0); u0=[0.85], ymax=ymax, color=blue_scale[1], la=0.8)
    plot_sdde_stable_bistability = plot_sdde_activation!(parameters_sdde_stable_bistability, t_pts; ymax=ymax, color=blue_scale[2], la=0.8)
    plot_sdde_single_pulse_response = plot_sdde_activation(parameters_sdde_single_pulse_response, t_pts, ymax=ymax)
    plot_sdde_intermediary_activation = plot_sdde_activation(parameters_sdde_intermediary_activation, t_pts, ymax=ymax)
    
    sdde_behaviour_recreation_figure = plot(plot_sdde_no_activation,plot_sdde_stochastic_pulsing,plot_sdde_oscillation,plot_sdde_stochastic_anti_pulsing,plot_sdde_heterogeneous_activation,plot_sdde_homogeneous_activation,plot_sdde_stochastic_switching,plot_sdde_stable_bistability,plot_sdde_single_pulse_response,plot_sdde_intermediary_activation,layout=(2,5),size=(3200,800))    
    make_figure(sdde_behaviour_recreation_figure, "SDDE_behaviour_recreation"; supp=true)
end

### Supplementary Figure 13-14, A1-A3 Correlation ###


# Shows correlation in a single simulation (Supplementary figure 13).
let 
    include("../Codebase/Functions_specific/A1_A3_correlation_analysis.jl")

    # Make simulation.
    parameters_sp = [9.,100.,10.0,0.05,2,0.05]
    sol = cle_activation(parameters_sp, (-1.0, 1000.0), 0.0, saveat=1.)
    
    p1 = plot(sol,idxs=[2,4],lw=4,color=[:orange :orangered4],label=["A1" "A3"],la=[0.9 0.7],ylimit=(0.0,0.7))
    s1 = getindex.(sol.u,2); s2 = getindex.(sol.u,4); 
    
    plot(xcorr(s1 .-mean(s1),s2 .-mean(s2)),xlimit=(0.0,2000),label="cross-correlation(A1,A3)",lw=7,la=0.7,ylimit=(-1,7.0)) 
    p2 = plot!(xcorr(s1 .-mean(s1),s1 .-mean(s1)),xlimit=(0.0,2000),label="cross-correlation(A1,A1)",lw=2,ylimit=(-1,7.0),color=:red,la=0.9)
    A1_A3_correlation_single_example_figures = plot(p1,p2,size=(1400,350),bottom_margin=7mm,layout = @layout [a b{0.4w}])
    
    make_figure(A1_A3_correlation_single_example_figures, "A1_A3_correlation_single_example"; supp=true)
end

# Shows correlation in a single simulation (Supplementary figure 14).
let 
    behaviours_to_params = get_behaviour_param_dict(dataset)
    sorted_list_of_behaviours = [:no_activation,:stochastic_pulsing,:oscillation,:stochastic_anti_pulsing,:homogeneous_activation,:heterogeneous_activation,:stochastic_switching,:stable_bistability,:single_response_pulse,:homogeneous_intermediate_activation]
    @time all_xcorr_plots = map(behaviour -> map(repeat -> plot_mean_xcorrs(filter(p -> 0.001<p[6]<0.1, behaviours_to_params[behaviour]), 100; l=50000.0), 1:5), sorted_list_of_behaviours);
    A1_A3_correlation_examples_figures = plot(vcat(all_xcorr_plots...)...,left_margin=-1mm,bottom_margin=-0.75mm,layout=(10,5),size=(1600,1600))
    
    make_figure(A1_A3_correlation_examples_figures, "A1_A3_correlation_examples"; supp=true)
end

### Supplementary Figure 15-16, Classifier Test ###

# Randomly selects 10 parameter sets from each behaviour.
# Need only be run once (to select the parameters). Does not need to be rerun when figures are regenerated.
# let
#     behaviours_to_params = get_behaviour_param_dict(dataset)
#     filter!(x -> 0.001 < x[2][6] < 0.1, behaviours_to_params)
#     randomly_selected_parameter_sets = Dict{Symbol,Vector{Vector{Float64}}}()
#     for behaviour in list_of_behaviours
#         randomly_selected_parameter_sets[behaviour] = rand(behaviours_to_params[behaviour], 10)
#     end
#     serialize("../Data/Simulated_data/$(dataset.dataset_tag)/randomly_selected_parameter_sets.jls", randomly_selected_parameter_sets)
# end

# General validation parameters.
begin 
    time_pts = (-1000.0, 0.0, 5000.0)
    ymax = 1.25; activation_lw=4;
    xguide = ""; yguide = "";
    idxs=[1]; lw=2; la=0.8; color=:darkblue;
    saveat=1.; maxiters=10000000;
    randomly_selected_parameter_sets = deserialize("../Data/Simulated_data/$(dataset.dataset_tag)/randomly_selected_parameter_sets.jls")
end

# No activation validation.
no_activation_validation_plot = let
    parameters_no_activation = randomly_selected_parameter_sets[:no_activation]
    initial_plots_no_activation = map(p -> plot_activation(p, time_pts; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_no_activation)
    plot(initial_plots_no_activation..., size=(1600,450), layout=(2,5))
end

# Stochastic pulsing validation.
stochastic_pulsing_validation_plot = let
    l_stochastic_pulsing = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]

    parameters_stochastic_pulsing = randomly_selected_parameter_sets[:stochastic_pulsing]
    initial_plots_stochastic_pulsing = map((p,l) -> plot_activation(p, (-l/5, 0.0, l); plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_stochastic_pulsing, l_stochastic_pulsing)
    plot(initial_plots_stochastic_pulsing..., size=(1600,450), layout=(2,5))
end

# Oscillation validation.
oscillation_validation_plot = let
    l_oscillation = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]

    parameters_oscillation = randomly_selected_parameter_sets[:oscillation]
    initial_plots_oscillation = map((p,l) -> plot_activation(p, (-l/5, 0.0, l); plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_oscillation, l_oscillation)
    plot(initial_plots_oscillation..., size=(1600,450), layout=(2,5))
end

# Stochastic anti pulsing validation.
stochastic_anti_pulsing_validation_plot = let
    l_stochastic_anti_pulsing = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]

    parameters_stochastic_anti_pulsing = randomly_selected_parameter_sets[:stochastic_anti_pulsing]
    initial_plots_stochastic_anti_pulsing = map((p,l) -> plot_activation(p, (-l/5, 0.0, l); plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_stochastic_anti_pulsing, l_stochastic_anti_pulsing)
    plot(initial_plots_stochastic_anti_pulsing..., size=(1600,450), layout=(2,5))
end

# Heterogeneous activation validation.
heterogeneous_activation_validation_plot = let
    l_heterogeneous_activation = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]
    lw_heterogeneous_activation = [2.5 ,2.5 ,2.5 ,2.5 ,1.5 ,2.5 ,2 ,2.5 ,2.5 ,2.5]

    parameters_heterogeneous_activation = randomly_selected_parameter_sets[:heterogeneous_activation]
    initial_plots_heterogeneous_activation = map((p, l, lw) -> plot_activations(p, (-l/5.0, 0.0, l), 3; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, lw=lw, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_heterogeneous_activation, l_heterogeneous_activation, lw_heterogeneous_activation)
    plot(initial_plots_heterogeneous_activation..., size=(1600,450), layout=(2,5))
end

# Homogeneous activation validation.
homogeneous_activation_validation_plot = let
    parameters_homogeneous_activation = randomly_selected_parameter_sets[:homogeneous_activation]
    initial_plots_homogeneous_activation = map(p -> plot_activation(p, time_pts; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_homogeneous_activation)
    plot(initial_plots_homogeneous_activation..., size=(1600,450), layout=(2,5))
end

# stochastic_switching validation.
stochastic_switching_validation_plot = let
    parameters_stochastic_switching = randomly_selected_parameter_sets[:stochastic_switching]
    initial_plots_stochastic_switching = map(p -> plot_activation(p, time_pts; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_stochastic_switching)
    plot(initial_plots_stochastic_switching..., size=(1600,450), layout=(2,5))
end

# Stable bistability validation.
stable_bistability_validation_plot = let
    function plot_stable_bistability_activation(p, l)
        plot_activations(p, (-l/5.0, 0.0, l), 1; plot_nc_minmaxes=true, colors=blue_scale[1], la=0.8, ymax=ymax, lw=3)
        plot_simulation!(p, (0.0, l), idxs=[1], color=blue_scale[3], ymax=ymax, la=0.5, u0=[1.0,1.0,1.0,1.0], lw=3, xticks=:auto, yticks=:auto)
        plot!(xlimit=(-l/5.0, l))
    end
    l_stable_bistability = [600., 600., 300., 600., 300., 300., 600., 300., 600., 300.]

    parameters_stable_bistability = randomly_selected_parameter_sets[:stable_bistability]
    initial_plots_stable_bistability = map((p,l) -> plot_stable_bistability_activation(p, l), parameters_stable_bistability, l_stable_bistability)
    plot(initial_plots_stable_bistability..., size=(1600,450), layout=(2,5))   
end

# Single response pulse validation.
single_response_pulse_validation_plot = let
    parameters_single_response_pulse = randomly_selected_parameter_sets[:single_response_pulse]
    initial_plots_single_response_pulse = map(p -> plot_activation(p, time_pts; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_single_response_pulse)
    plot(initial_plots_single_response_pulse..., size=(1600,450), layout=(2,5))
end

# Intermediary activation validation.
intermediary_activation_validation_plot = let
    parameters_intermediary_activation = randomly_selected_parameter_sets[:homogeneous_intermediate_activation]
    initial_plots_intermediary_activation = map(p -> plot_activation(p, time_pts; plot_nc_minmaxes=true, xguide=xguide, yguide=yguide, xticks=:auto, yticks=:auto, idxs=idxs, lw=lw, la=la, colors=color, ymax=ymax, saveat=saveat, activation_lw=activation_lw, maxiters=maxiters), parameters_intermediary_activation)
    plot(initial_plots_intermediary_activation..., size=(1600,450), layout=(2,5))
end

# Saves the example classification figures (split into two).
let 
    example_trajectories_plot_1 = plot(no_activation_validation_plot, stochastic_pulsing_validation_plot, oscillation_validation_plot, stochastic_anti_pulsing_validation_plot, heterogeneous_activation_validation_plot,
                                 left_margin=-1mm,bottom_margin=-0.75mm,layout=(5,1),size=(1600,1600))
    make_figure(example_trajectories_plot_1,"Example_classifications"; tag="1", supp=true)
end
let 
    example_trajectories_plot_2 = plot(homogeneous_activation_validation_plot, stochastic_switching_validation_plot, stable_bistability_validation_plot, single_response_pulse_validation_plot, intermediary_activation_validation_plot,
                                    left_margin=-1mm,bottom_margin=-0.75mm,layout=(5,1),size=(1600,1600))
    make_figure(example_trajectories_plot_2,"Example_classifications"; tag="2", supp=true)
end


### Supplementary Figure 17, Active/Inactive State Definition ###
let 
    # Nullcline variation with S.
    xmax = 1.; ymax=2.;
    parameters = [2.25,1.1,0.1,0.05,2,0.05]
    ncmm = nc_local_min_max(parameters)
    plot_nullcline_sets(parameters,:S,2.0:0.20:2.8,colors=make_colour_scale(RGB{Float64}(0.2,0.,0.),RGB{Float64}(1.0,0.,0.),5),xmax=xmax,ymax=ymax,legendfontsize=12)
    p1 = plot!([[ncmm[1],ncmm[1]] [ncmm[2],ncmm[2]]],[0.,ymax],xticks=:auto,yticks=:auto,color=:black,lw=3,la=0.7,linestyle=:dot,label="")

    # Defining states in simulation.
    parameters = [2.25,1.1,10.0,0.05,2,0.05];
    ncmm = nc_local_min_max(parameters)
    sT = 200.
    leng = 1000.
    ymax = 1.25
    ia_colors = [RGB{Float64}(0.,0.7,0.7) RGB{Float64}(0.,.0,1.) RGB{Float64}(0.7,0.0,0.7)];
    function find_ia_color(val,ncmm1,ncmm2)
        (val[1]<ncmm1) && return ia_colors[1]
        (val[1]>ncmm2) && return ia_colors[3]
        return ia_colors[2]
    end

    sol = cle_activation(parameters, (0.0,1000.0), 200.0, saveat=0.1)
    lPre = Int64(sT/0.1); lPost = length(sol.u)-Int64(sT/0.1);
    colors = find_ia_color.(sol.u,ncmm...);

    plot_nc_min_maxes(args...;kwargs...) = (plot(); plot_nc_min_maxes!(args...;kwargs...);)
    function plot_nc_min_maxes!(p_vals::Vector{Float64}, tspan::Tuple{Float64,Float64}; model=model3, lw=1, la=0.9, plot_label=false)
        ncmm = nc_local_min_max(p_vals); 
        (length(ncmm) == 0) && push!(ncmm,large_v0_inactive_active_threshold(p_vals[5]))
        foreach(val -> plot!([tspan...],[val,val],color=:black,linestyle=:dot,lw=lw,la=la,label=""), ncmm)
        plot_label && plot!([],[]; color=:black, linestyle=:dot, lw=lw, la=la, label="Inactive/active state thresholds")
        plot!()
    end

    plot(sol,idxs=[1],label="",lw=4,la=0.6,color=colors)
    plot_nc_min_maxes!(parameters, (0.0,1000.0); plot_label=true)
    plot_vertical_line!(sT,(ymax==Inf ? maximum(first.(sol.u)) : max(ymax,maximum(first.(sol.u)))),lw=3, label="System activation")
    p2 = plot!([0.,0.],[[0.,0.] [0.,0.] [0.,0.]],color=ia_colors,lw=4,la=0.6,label=["Inactive state" "Intermediary state" "Active state"],xguide="Time",yguide="σ",legend=:topright,legendfontsize=9,xlimit=(0.,leng),ylimit=(0.,ymax),xticks=[],yticks=[])

    # Defining states in phase space.
    plot(sol,idxs=(1,4),label="",lw=2,la=0.4,color=colors)
    plot!([[ncmm[1],ncmm[1]] [ncmm[2],ncmm[2]]],[0.,ymax],color=:black,lw=3,la=0.7,linestyle=:dot,label=["Inactive/active state thresholds" ""])
    plot_nullcline_set!(parameters,ymax=ymax,xticks=[],yticks=[],left_margin=4mm)
    p3 = plot!([0.,0.],[[0.,0.] [0.,0.] [0.,0.]],color=ia_colors,lw=4,la=0.6,label=["Inactive state" "Intermediary state" "Active state"],xguide="σ",yguide="A",legend=:topright,legendfontsize=9,xlimit=(0.,ymax),ylimit=(0.,ymax),xticks=[],yticks=[])

    # Plots various nullclines for various v0.
    function plot_ncmm!(parameters,ymax;color=:black,lw=3,linestyle=:dot,la=0.7,label="")
        ncmm = nc_local_min_max(parameters)
        if length(ncmm) == 2
            return plot!([[ncmm[1],ncmm[1]] [ncmm[2],ncmm[2]]],[0.,ymax],color=color,lw=lw,la=la,linestyle=linestyle,label=[label ""])
        end
        plot!([large_v0_inactive_active_threshold(2),large_v0_inactive_active_threshold(2)],[0.,ymax],color=color,lw=lw,la=la,linestyle=linestyle,label=label)
    end
    ymax=5

    colors = make_colour_scale(RGB{Float64}(1.,0.,0.), RGB{Float64}(0.,0.,1.), 6);
    v0s = [0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.14,0.15]
    parameters = [1.65,.2,0.1,0.08,2,0.]

    plot();
    for i = 1:6
        plot_nullcline_set!(updated_p(parameters, :v0, v0s[i]), ymax=ymax, xticks=:auto, yticks=:auto, left_margin=4mm, color2=colors[i])
        plot_ncmm!(updated_p(parameters, :v0, v0s[i]), ymax; color=colors[i])
    end
    for i = 7:9
        plot_nullcline_set!(updated_p(parameters, :v0, v0s[i]),ymax=ymax,xticks=[],yticks=[],left_margin=4mm,color2=colors[6])
    end
    p4 = plot!([],fill([],length(colors)),xlimit=(0.,0.85),ylimit=(0.,ymax),color=hcat(colors...),label=["v0=$(v0s[1])" "v0=$(v0s[2])" "v0=$(v0s[3])" "v0=$(v0s[4])" "v0=$(v0s[5])" "v0=$(v0s[6]),$(v0s[7]),$(v0s[8]),$(v0s[9])"])

    definting_active_inactive_plot = plot(p1, p2, p4, p3, size=(1200,700), left_margin=12mm, bottom_margin=7mm, layout=(2,2))
    make_figure(definting_active_inactive_plot,"Definting_active_inactive"; supp=true)
end


### Supplementary Figure 18-29, Line Parameter Dependency Analysis ###

# Finds and saves all transition points.
# let 
#     include("../Codebase/Functions_specific/transition_line_analysis.jl")
#     transition_points = find_transition_points(dataset)
#     serialize("Data/Simulated_data/$(dataset.dataset_tag)/transition_points.jls",transition_points);    
# end

# Prints data of all line transitions.
let 
    include("../Codebase/Functions_specific/transition_line_analysis.jl")
    transition_points = deserialize("../Data/Simulated_data/$(dataset.dataset_tag)/transition_points.jls")
    sorted_transitions, all_linedata = get_sorted_transitions(transition_points)
    display_transitions(sorted_transitions, all_linedata)    
end

# Prints the fitted lines on behaviour map.
let
    gr()
    params = [5.0,0.05,3,0.1];
    bg = BehaviourGrid(params,dataset);
    behaviour_grid_plot = plot!(plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="",yguide=""),size=(600,400))

    include("../Codebase/Functions_specific/transition_line_analysis.jl")
    transition_points = deserialize("../Data/Simulated_data/$(dataset.dataset_tag)/transition_points.jls")
    sorted_transitions, all_linedata = get_sorted_transitions(transition_points)
    for idx = 1:7
        trans = sorted_transitions[idx][1]
        line = all_linedata[trans][findfirst(dataset.τ_grid.==5.0),findfirst(dataset.v0_grid.==0.05),findfirst(dataset.n_grid.==3),findfirst(dataset.η_grid.==0.1)]
        tps = transition_points[trans][findfirst(dataset.τ_grid.==5.0),findfirst(dataset.v0_grid.==0.05),findfirst(dataset.n_grid.==3),findfirst(dataset.η_grid.==0.1)]
        plot!(x->line.a*x+line.b,minimum(first.(tps)),maximum(first.(tps)),color=(idx+2),la=0.8,xlimit=(0.1,100.0),ylimit=(1.0,100.0),background_color=:transparent,subplot_background_color=:transparent,xaxis=:log10,yaxis=:log10,xguide="",yguide="",xticks=[],yticks=[],size=(600,400),label="",lw=5, inset = (1, bbox(0.,0.0,1.0,1.0)), subplot = 2)
        println("\n",trans)
        println(line.a,"          ",line.b)
    end
    transition_lines_on_behaviour_map_figure = plot!(xticks=[],yticks=[])
    make_figure(transition_lines_on_behaviour_map_figure,"Transition_lines_on_behaviour_map"; supp=true)
end

# Prints the seven most common transitions.
let
    include("../Codebase/Functions_specific/transition_line_analysis.jl")
    transition_points = deserialize("../Data/Simulated_data/$(dataset.dataset_tag)/transition_points.jls")
    sorted_transitions, all_linedata = get_sorted_transitions(transition_points)
    for transition in sorted_transitions[1:7]
        parameter_line_plot = plot_lines_all_params(transition[1], dataset, all_linedata, transition_points)
        make_figure(parameter_line_plot,"Transition_Line_Parameter_Dependence"; supp=true, tag="$(transition[1])___$(transition[2])")
    end    
end

# Plots the transition lines in linear space.
let 
    gr()
    dataset_lin = DataSet("LinearMap",range(0.1,stop=10.0,length=1000),range(0.1,stop=10.0,length=1000),[5.0],[0.05],[3.0],[0.1])
    bg_lin = BehaviourGrid([5.0,0.05,3.0,0.1], dataset_lin)
    foreach(i -> (bg_lin.behaviours[371,i]==:heterogeneous_activation) && (bg_lin.behaviours[371,i]=:homogeneous_activation), 1:1000)
    lin_behaviour_grid_plot = plot!(plot_behaviour_grid(bg_lin,start_s_slice=100,xguide="",yguide="",xticks=[],yticks=[]),size=(1200,800))

    DSs = Vector{Tuple{Float64,Float64}}()
    for (Di,D) in enumerate(dataset_lin.D_grid[1:end]), (Si,S) in enumerate(dataset_lin.S_grid[1:end-1])
        if (bg_lin.behaviours[Si,Di]==:no_activation) && (bg_lin.behaviours[Si+1,Di]==:stable_bistability)
        push!(DSs,(D,S)) 
        end
    end
    xi = first.(DSs); yi = last.(DSs);    
    A = [xi ones(length(xi))]
    a, b = A\yi
    er = sum((yi .- (a*xi .+ b) ) .^2)/length(DSs);
    println(a,"\t", b, "\t", er)

    scatter(DSs,color=RGB{Float64}(0,0.25,0),ms=5)
    NA_SB_line_plot = plot!(x->a*x+b,0.07,0.9,color=1,lw=8,la=0.5,xlimit=(0.08,0.75),ylimit=(2.08,2.44),legend=:none,size=(600,200),xticks=[0.1,0.4,0.7],yticks=[2.1,2.2,2.3,2.4],xguide="D", yguide="S")

    behaviour_map_linear_space_plot = plot(lin_behaviour_grid_plot,p0_gr ,NA_SB_line_plot, p0_gr, size=(1200,1000),layout=@layout[a{0.80h}; b c{0.60w} d])
    make_figure(behaviour_map_linear_space_plot,"Behaviour_map_linear_space"; supp=true)
end


### Supplementary Figure 30, Bifurcation of Bifrucation Points ###


### Unsorted Supplementary Figures ###

# Generates example simulations of initial 4 behaviour types.
let
    t_pts = (-100.0, 0.0, 500.0)
    ymax = 1.25;

    # Parameters.
    parameters_single_response_pulse =      [2.2,  6.0,  10.0, 0.1,   2, 0.025]
    parameters_stochastic_pulsing =         [2.8,  10.0, 10.0, 0.05,  2, 0.025]
    parameters_heterogeneous_activation =   [2.25, 0.1,  10.0, 0.05,  2, 0.025]
    parameters_homogeneous_activation =     [5.0,  0.1,  10.0, 0.05,  2, 0.025]

    # Plots.
    plot_homogeneous_activation = plot_activation(parameters_homogeneous_activation, t_pts; idxs = [1], ymax=ymax)
    plot_single_response_pulse = plot_activation(parameters_single_response_pulse, t_pts; idxs = [1], ymax=ymax)
    plot_heterogeneous_activation = plot_activations(parameters_heterogeneous_activation, t_pts, 3; ymax=ymax)
    plot_stochastic_pulsing = plot_activation(parameters_stochastic_pulsing, t_pts; idxs = [1], ymax=ymax)

    # Make plot.
    four_behaviour_example_sim_plot = plot(plot_homogeneous_activation,plot_single_response_pulse,plot_heterogeneous_activation,plot_stochastic_pulsing,layout=(2,2),size=(1400,800))    
    make_figure(four_behaviour_example_sim_plot, "Four_behaviour_example_sims"; supp=true)
end