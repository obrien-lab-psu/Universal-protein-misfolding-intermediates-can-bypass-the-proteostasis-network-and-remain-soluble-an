package require topotools
display rendermode GLSL
axes location off

color Display {Background} white

mol new ./state_3.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 10
mol selection {all}
mol material AOChalky
mol addrep top
