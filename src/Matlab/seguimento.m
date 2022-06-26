function rumo_final=seguimento(in)

global j;
global Path;
global rumo;
global vec_anterior;
global Influencia;

len=length(Path);

while j<len
    if norm([in(2) in(1)]-[Path(j,1) Path(j,2)])>Influencia-1
        break
    end
  j=j+1;
end

if j<len

    vec_dest=[Path(j,1) Path(j,2)];
    vec_actual=[in(2) in(1)];

    dir_dest=vec_dest-vec_actual;
    dir_dest=dir_dest/norm(dir_dest);

    dir_actual=vec_actual-vec_anterior;
    dir_actual=dir_actual/norm(dir_actual);

    rota=acos(dir_actual*dir_dest');
    vertical=acos(dir_actual*[0; 1]);

    sentido_vertical=det(cat(1,dir_actual,[0 1]));
    sentido_rot=det(cat(1,dir_actual,dir_dest));

    if sentido_vertical<0
        vertical=-vertical;
    end

    if sentido_rot>0
        rota=-rota;
    end

    rumo=rota+vertical;
    vec_anterior=vec_actual;
end

disp(rumo)
if isnan(rumo)
    rumo = 0;
end
rumo_final=rumo;
