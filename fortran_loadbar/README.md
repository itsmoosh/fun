Prints a loading bar to the terminal that is self-overwriting.

Output cursor is stepped backward after each time the loading
bar is updated, unless percent >= 100. Call loadbar(100) after
the last sub-100 update in order to place the cursor at a normal
location.

Inputs:
	percent -- integer in the range [0, 100]

Usage:
	integer :: load_pct = 95
	call loadbar(load_pct)
