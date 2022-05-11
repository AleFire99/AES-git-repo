within AES_project_2021_2022;

function NotteGiorno
input Integer t;  
output Integer R;

protected 

Real hour:=t/3600;
Real gg:=hour/24;
Real hourloc:=hour-gg*24;

algorithm
    if (hourloc<=8 and hourloc>=22) then R:=0; end if;
end NotteGiorno;
