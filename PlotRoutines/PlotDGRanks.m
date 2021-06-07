% script for plotting the ranks of blocks on a DG stiffness matrix
% written by Boris Bonev 12/2018
clear all
close all

%% choose parameters
h = 4;
p = 2;
rktol = 1e-6;
rkmode = 'relative rank';

%% generate 2D matrix
GenMatrixPoissonSquare2D('test.mat', h, p);
load('test.mat');

% get dimensions 2D
nCells = h^2*2;
nPoints = size(A,2)/nCells;

%% generate 3D matrix
% GenMatrixPoissonSquare3D('test.mat', h, p);
% load('test.mat');
% 
% % dimensions in 3D
% nCells = h^3*5;
% nPoints = size(A,2)/nCells;

%% get DG blocks as separators
% separate into individual cells
seps = {};
for itet = 1:nCells
    seps{1,itet} = (1:nPoints)+(itet-1)*nPoints;
end
% get inverse
B = inv(A);

%% plot cartesian product blocks
% % plot without the reordering
% figure()
% spy(A)
% PlotRanksReorder(A,seps,rkmode, rktol)
% 
% % plot after reordering
% PlotSparsityReorder(A, sep_tree, 5);
% PlotRanksReorder(A,{sep_tree{1,:}},rkmode, rktol)
%  
% % plot ranks of the inverse
% figure()
% spy(B)
% PlotRanksReorder(B,seps,rkmode, rktol)
% 
% % plot B after reordering
% PlotSparsityReorder(B, sep_tree, 5);
% PlotRanksReorder(B,{sep_tree{1,:}},rkmode, rktol)

%% plot block cluster ranks
%seps = {[seps{1:16}], [seps{16+1:end}]};
cluster_tree = GenClusterTreeBisection(seps);
PlotRanksClusterTree(A, cluster_tree, rkmode, rktol)
caxis([0 1])
PlotRanksClusterTree(B, cluster_tree, rkmode, rktol)
caxis([0 1])

[cluster_tree, perm] = Sep2ClusterTree(sep_tree);
PlotRanksClusterTree(A(perm, perm), cluster_tree, rkmode, rktol)
caxis([0 1])
PlotRanksClusterTree(B(perm, perm), cluster_tree, rkmode, rktol)
caxis([0 1])