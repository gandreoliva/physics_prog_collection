# Structure of a static star

Qualitative solution to the stellar structure equations with the shooting method. No stellar evolution.

## Usage

Compile with `make`, example input text files are in the folder `stars` and contain the stellar mass (solar masses), luminosity (solar luminosities), temperature (K), hydrogen fraction (ratio to solar), metallicity (ratio to solar metallicity).

Run with, e.g. `./static_star.bin  < stars/sun.txt`.