% Wrapper only: the author's three SDP functions are called without edits.
root=getenv('EXTERNAL_WORKROOT'); out=getenv('EXTERNAL_NATIVE_OUT');
addpath(genpath([root '/matlab_runtime/YALMIP_release']));
addpath([root '/matlab_runtime/sdpt3_release']);
addpath([root '/matlab_runtime/sdpt3_release/HSDSolver']);
addpath(genpath([root '/matlab_runtime/sdpt3_release/Solver']));
repo=[root '/author_code/RDM_Constraints_Renormalization'];
addpath([repo '/auxFunctions']);addpath([repo '/spinChains_1D/MPS_Relaxation']);
load([out '/input.mat'],'H','V','L','R','last_k');
k0=3;D=2;CGmaps=struct('V0',V,'L',{cell(last_k,1)},'R',{cell(last_k,1)});
for k=k0:last_k; CGmaps.L{k}=L;CGmaps.R{k}=R;end
ops=sdpsettings('solver','sdpt3','verbose',1,'debug',1,'sdpt3.gaptol',1e-10,'sdpt3.inftol',1e-10,'sdpt3.maxit',150);
ops.sdpt3.model=1; % Explicit standard solver; automatic HSD failed the probe.
t=tic;
[energy,vars,info]=Relaxation_dual_SDP(last_k,D,k0,H,CGmaps,ops);
elapsed=toc(t);
duals={value(vars.alpha),value(vars.betaL),value(vars.betaR)};
for k=k0:last_k-1;duals{end+1}=value(vars.gammaL{k});duals{end+1}=value(vars.gammaR{k});end
r=struct('matlab',version,'yalmip',yalmip('version'),'author_function',which('Relaxation_dual_SDP'),...
 'energy',energy,'duals',{duals},'candidate_seconds',elapsed,'solver_seconds',info.solvertime,...
 'model_seconds',info.yalmiptime,'problem',info.problem,'info',info.info);
fid=fopen([out '/matlab_output.json'],'w');fprintf(fid,'%s',jsonencode(r,PrettyPrint=true));fclose(fid);
