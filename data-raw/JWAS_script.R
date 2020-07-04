# the script
JWAS_script <-
"using JWAS,CSV,DataFrames,DelimitedFiles;

phenofile  = 'ie.txt';
genofile   = 'ig.txt';

phenotypes = CSV.read(phenofile,delim = ' ',header=true);

model_equation1  = 'phenotypes = intercept';

R      = Meta.parse(ARGS[3]); #residual variance
model1 = build_model(model_equation1,R);
G      = Meta.parse(ARGS[4]); #genetic variance

add_genotypes(model1,genofile,G,separator=' ');

brn = Meta.parse(ARGS[2]);
chain = Meta.parse(ARGS[1]);

println(R);
println(G);

outputMCMCsamples(model1,'phenotype')
outputMCMCsamples(model1,'genetic_variance')


if ARGS[6] != 'false'
    out1=runMCMC(model1,phenotypes,methods=ARGS[5],chain_length=chain,burnin=brn,printout_frequency=10000,output_samples_frequency=50,Pi=Meta.parse(ARGS[6]));
else
    out1=runMCMC(model1,phenotypes,methods=ARGS[5],chain_length=chain,burnin=brn,printout_frequency=10000,output_samples_frequency=50,estimatePi=true);
end


# get h and write results

println('MCMC finished. Estimating h.');
res_gv=JWAS.misc.get_additive_genetic_variances(model1, 'MCMC_samples_marker_effects_phenotypes.txt');
res_rv=open(readdlm,'MCMC_samples_residual_variance.txt');
res_rv=res_rv[2:length(res_rv),1];
res_rv=convert(Array{Float64,1},res_rv);
h=JWAS.misc.get_heritability(res_gv, res_rv);

writedlm('est_effects.txt', out1['Posterior mean of marker effects']);
writedlm('h.txt',mean(h));"

# fix single -> double quotes
analysis.jl <- gsub("\'", "\"", JWAS_script)

# save
usethis::use_data(analysis.jl, overwrite = TRUE, internal = T)


