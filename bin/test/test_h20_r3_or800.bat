cd..
BeamSplitting -p 10 100 40 20 -ri 1.31 -b 0 800 -g 0 801 -t 180 -rn 3 -o %~n0


goto endcomment

######################################

*** Command line arguments ***

Note: all of these arguments must be present except ones marked as optional.


List of argumentst:

-p type hh rad (cd) - particle params: 

	type - type of particle (int) (0 - Hexagonal prism, 10 - Concave hexagonal prism);
	hh - half height (double);
	rad - radius (double);
	cd - cavity depth (double) (if type of particle is '10');


-ri re - real part of refraction index (double);


-rn num - number of the reflections (int);


-b beg end - beta number of orientations:

	beg - start orientation (int);
	end - end orientation (int);


-g beg end - gamma number of orientations:

	beg - start orientation (int);
	end - end orientation (int) (must be odd);


-t size - other tracing params:
	size - size of the output array (int);


-r - random orientation of particle due tracing (otherwise orientation is fixed) (optional);


-o filename - output file name (string) (optional);

######################################

:endcomment