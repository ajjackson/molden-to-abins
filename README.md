# Molden to Abins format converter

This script is intended to facilitate simulation of inelastic neutron
scattering spectra from a file in Molden format.


## Usage

By default, the script prints to standard output. From a unix-like environment this can be redirected to make a new .json file

```
molden_to_abins.py my_input.mol > my_output.json
```

or alternatively the ``--output`` parameter can be used.


```
molden_to_abins.py my_input.mol --output my_output.json
```

The resulting JSON file can be used with Mantid by selecting "JSON" as
the "AbInitioProgram" in the Abins or Abins2D algorithm. This requires
Mantid version 6.10 or later, or a Nightly build from April/May 2024.

## Requirements

The script was developed with Python 3.10 and requires
[ASE](https://wiki.fysik.dtu.dk/ase/); this is not
included with Mantid so you may need to install it with `pip install
ase`.

## Limitations

- The correctness has not been fully verified yet against another
  calculation method. It _appears_ that the eigenvectors are
  normalised correctly but validation is recommended before use in
  production.

- The Molden format does not include atomic masses, so we get the
  standard masses for each element (i.e. standard isotopic mixtures)
  using the ASE library.
