function VisualizeTreeCG2D(elim_tree)

Globals2D;
GlobalsCG2D;

v = zeros(nTotal,1);
for i=1:size(elim_tree,2)
  v(elim_tree{1,i}) = i;
end

SurfCG2D(v)

end