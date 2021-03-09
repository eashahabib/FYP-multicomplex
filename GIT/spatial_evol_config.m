function gp = spatial_evol_config(gp)


gp.fitness.fitfun = @quartic_fitfun;

%run control
gp.runcontrol.pop_size = 250;
gp.runcontrol.runs = 1;
gp.runcontrol.timeout = 60;
gp.runcontrol.parallel.auto = true; %parallel computing if available
gp.state.force_compute_theta = false;

%selection
gp.selection.tournament.size = 20;
gp.selection.tournament.p_pareto = 0.3;
gp.selection.elite_fraction = 0.3;

%genes
gp.genes.max_genes = 10;
gp.genes.multigene = false;

%training grid
[x1, x2] = meshgrid(-5:0.4:5, -5:0.4:5);
x1 = x1(:); x2 = x2(:);
y = 1./(1+x1.^-4) + 1./(1+x2.^-4);
gp.userdata.ytrain = y;
gp.userdata.xtrain = [x1 x2];

%test grid
[x1, x2] = meshgrid(-5:0.2:5, -5:0.2:5);
x1 = x1(:); x2 = x2(:);
y = 1./(1+x1.^-4) + 1./(1+x2.^-4);
gp.userdata.ytest = y;
gp.userdata.xtest = [x1 x2];

%function nodes
gp.nodes.functions.name = {'times' 'minus' 'plus' 'rdivide' 'square' 'sin' 'cos' 'exp' 'mult3' 'add3' 'sqrt' 'cube' 'power' 'negexp' 'neg' 'abs' 'log'};
