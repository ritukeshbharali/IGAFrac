function writeVtk(npoin,nelem,nnode,lnods,coord,istep,cont1,displacement,folder_name)

fname = [folder_name,'/time',num2str(istep),'.vtk'];
out = fopen(fname,'w');

% Start writing:
% Header
fprintf(out,'# vtk DataFile Version 2.0\n');
fprintf(out,'time_10.vtk\n');
fprintf(out,'ASCII\n');
fprintf(out,'DATASET UNSTRUCTURED_GRID\n');

% Write nodal coordinates:
fprintf(out,'POINTS %5d float\n',npoin);

dummy = 0.0;
for ipoin=1:npoin
    fprintf(out,'%14.6f %14.6f %14.6f\n',coord(ipoin,1),coord(ipoin,2),dummy);
end

% Write element connectivity:
iconst1 = nelem*(nnode+1);

fprintf(out,'CELLS %5d %5d\n', nelem,iconst1);
for ielem=1:nelem
    fprintf(out,'%5d',nnode);
    for inode=1:nnode
        fprintf(out,' %5d',(lnods(ielem,inode)-1));
    end
    fprintf(out,'\n');
end

% Write cell types:
if(nnode == 8)
    ntype = 23;
end

if(nnode == 4)
    ntype = 9;
end

if(nnode == 3)
    ntype = 5;
end

fprintf(out,'CELL_TYPES %5d\n',nelem);
for i=1:nelem
    fprintf(out,'%2d\n', ntype);
end

% Write the displacement as scalars
fprintf(out,'POINT_DATA %5d\n',npoin);
fprintf(out,'VECTORS U float \n');

%fprintf(out,'LOOKUP_TABLE default\n');
for ipoin = 1:npoin
    fprintf(out,'%14.6f %14.6f %14.6f\n',displacement(ipoin,1),displacement(ipoin,2),dummy);
end

fprintf(out,'SCALARS Con  float  1\n');
fprintf(out,'LOOKUP_TABLE default\n');
for ipoin=1:npoin
    fprintf(out,'%14.6e\n',cont1(ipoin));
end
fclose(out);

end %endfunction
