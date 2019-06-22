function impatch(V)

figure; 
data = smooth3(V,'box',5);
patch(isocaps(data,.5),...
'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data,.5),...
'FaceColor','blue','EdgeColor','none');
isonormals(data,p1)
view(3);
axis vis3d tight
camlight left;
colormap jet
lighting gouraud
end