J/MNRAS/423/122   Abundances of 93 solar-type Kepler targets    (Bruntt+, 2012)
================================================================================
Accurate fundamental parameters and detailed abundance patterns from
spectroscopy of 93 solar-type Kepler targets.
    Bruntt H., Basu S., Smalley B., Chaplin W.J., Verner G.A., Bedding T.R.,
    Catala C., Gazzano J.-C., Molenda-Zakowicz J., Thygesen A.O.,
    Uytterhoeven K., Hekker S., Huber D., Karoff C., Mathur S., Mosser B.,
    Appourchaux T., Campante T.L., Elsworth Y., Garcia R.A., Handberg R.,
    Metcalfe T.S., Quirion P.-O., Regulo C., Roxburgh I.W., Stello D.,
    Christensen-Dalsgaard J., Kawaler S.D., Kjeldsen H., Morris R.L.,
    Quintana E.V., Sanderfer D.T.
   <Mon. Not. R. Astron. Soc., 423, 122-131 (2012)>
   =2012MNRAS.423..122B
================================================================================
ADC_Keywords: Stars, dwarfs ; Stars, G-type ; Abundances
Keywords: stars: abundances - stars: atmospheres -
          stars: fundamental parameters - stars: solar-type

Abstract:
    We present a detailed spectroscopic study of 93 solar-type stars that
    are targets of the NASA/Kepler mission and provide detailed chemical
    composition of each target. We find that the overall metallicity is
    well represented by Fe lines. Relative abundances of light elements
    (CNO) and {alpha} elements are generally higher for low-metallicity
    stars. Our spectroscopic analysis benefits from the accurately
    measured surface gravity from the asteroseismic analysis of the Kepler
    light curves. The accuracy on the log g parameter is better than
    0.03dex and is held fixed in the analysis.

Description:
    The spectra were obtained with the ESPaDOnS spectrograph at the 3.6-m
    Canada-France-Hawaii Telescope (CFHT) in USA and with the NARVAL
    spectrograph mounted on the 2-m Bernard Lyot Telescope at the Pic du
    Midi Observatory in France. In both the facilities, the observations
    were carried out as service observations from May to September in
    2010.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table1.dat     37     1042   Atomic data used in the spectral analysis
table3.dat     96       93   Observed targets and their properties
table4.dat    218      186   Abundances relative to the Sun and the number of
                             lines used for 13 elements
--------------------------------------------------------------------------------

See also:
   V/133 : Kepler Input Catalog (Kepler Mission Team, 2009)

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   2-  5  A4    ---     El        Element with ionization state
   8- 15  F8.3  0.1nm   lambda    wavelength {lambda} ({AA})
  17- 22  F6.3  eV      Ep        excitation potential
  24- 29  F6.3  [-]     loggfV    Oscillator strength from VALD
  31- 37  F7.3  [-]     loggf     Adjusted oscillator strength value
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  8  I8    ---     KIC     KIC number (Cat. V/133)
  10- 14  I5    ---     HIP     ? HIP number (Cat. I/311)
  16- 21  I6    ---     HD      ? HD number (Cat. III/135)
  23- 28  F6.3  mag     VTmag   Tycho V magnitude
  30- 34  F5.3  mag     VT-K    VT-Ks colour index
  36- 39  F4.2  mag     E(B-V)  ? B-V colour excess
  41- 44  F4.2  mag   e_E(B-V)  ? rms uncertainty on E(B-V)
  46- 49  I4    K       Teff1   ? Effective temperature  from KIC
  51- 54  F4.2  [cm/s2] logg1   ? Gravity surface from KIC
  56- 60  F5.2  [Sun]   [Fe/H]1 ? Metallicity from KIC
  62- 65  I4    K       Teff2   Effective temperature from spectroscopy
                                 (+/-70K)
  67- 70  F4.2  [cm/s2] logg2   Gravity surface from spectroscopy (+/-0.08dex)
  72- 75  I4    K       Teff    Final (asteroseismic) effective temperature
                                 (+/-60K)
  77- 80  F4.2  [cm/s2] logg    Final (asteroseismic) gravity surface
                                 (+/-0.03dex)
  82- 86  F5.2  [-]     [Fe/H]  Final (asteroseismic) metallicity (+/-0.06dex)
  88- 91  F4.2  km/s    xi      Final (asteroseismic) microturbulent velocity
                                 {xi} (+/-0.06km/s) (1)
  93- 96  F4.1  km/s    vsini   Final (asteroseismic) rotational velocity (1)
--------------------------------------------------------------------------------
Note (1): We only list the `asteroseismic' values of vsini and xi, since
          the `spectroscopic' values are almost identical.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table4.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label  Explanations
--------------------------------------------------------------------------------
   1-  8  I8    ---     KIC    KIC number (1)
      10  I1    ---     Ion    [1/2] Ionisation state:
                               1 for neutral, 2 for singly-ionized
  12- 16  F5.2  [Sun]   Li     ?=- Lithum abundance (2)
  19- 22  F4.2  [Sun] e_Li     ?=- rms uncertainty on Li
  24- 26  I3    ---   o_Li     ?=- Number of lines used for Li
  28- 32  F5.2  [Sun]   C      ?=- Carbon abundance (2)
  35- 38  F4.2  [Sun] e_C      ?=- rms uncertainty on C
  40- 42  I3    ---   o_C      ?=- Number of lines used for C
  44- 48  F5.2  [Sun]   N      ?=- Nitrogen abundance (2)
  51- 54  F4.2  [Sun] e_N      ?=- rms uncertainty on N
  56- 58  I3    ---   o_N      ?=- Number of lines used for N
  60- 64  F5.2  [Sun]   O      ?=- Oxygen abundance (2)
  67- 70  F4.2  [Sun] e_O      ?=- rms uncertainty on O
  72- 74  I3    ---   o_O      ?=- Number of lines used for O
  76- 80  F5.2  [Sun]   Na     ?=- Sodium abundance (2)
  83- 86  F4.2  [Sun] e_Na     ?=- rms uncertainty on Na
  88- 90  I3    ---   o_Na     ?=- Number of lines used for Na
  92- 96  F5.2  [Sun]   Mg     ?=- Magnesium abundance (2)
  99-102  F4.2  [Sun] e_Mg     ?=- rms uncertainty on Mg
 104-106  I3    ---   o_Mg     ?=- Number of lines used for Mg
 108-112  F5.2  [Sun]   Si     ?=- Silicium abundance (2)
 115-118  F4.2  [Sun] e_Si     ?=- rms uncertainty on Si
 120-122  I3    ---   o_Si     ?=- Number of lines used for Si
 124-128  F5.2  [Sun]   Ca     ?=- Calcium abundance (2)
 131-134  F4.2  [Sun] e_Ca     ?=- rms uncertainty on Ca
 136-138  I3    ---   o_Ca     ?=- Number of lines used for Ca
 140-144  F5.2  [Sun]   Ti     ?=- Titanium abundance (2)
 147-150  F4.2  [Sun] e_Ti     ?=- rms uncertainty on Ti
 152-154  I3    ---   o_Ti     ?=- Number of lines used for Ti
 156-160  F5.2  [Sun]   V      ?=- Vanadium abundance (2)
 163-166  F4.2  [Sun] e_V      ?=- rms uncertainty on V
 168-170  I3    ---   o_V      ?=- Number of lines used for V
 172-176  F5.2  [Sun]   Cr     ?=- Chromium abundance (2)
 179-182  F4.2  [Sun] e_Cr     ?=- rms uncertainty on Cr
 184-186  I3    ---   o_Cr     ?=- Number of lines used for Cr
 188-192  F5.2  [Sun]   Fe     ?=- Iron abundance (2)
 195-198  F4.2  [Sun] e_Fe     ?=- rms uncertainty on Fe
 200-202  I3    ---   o_Fe     ?=- Number of lines used for Fe
 204-208  F5.2  [Sun]   Ni     ?=- Nickel abundance (2)
 211-214  F4.2  [Sun] e_Ni     ?=- rms uncertainty on Ni
 216-218  I3    ---   o_Ni     ?=- Number of lines used for Ni
--------------------------------------------------------------------------------
Note (1): For each star there are two rows: One for the neutral species and
          one for the singly-ionized species (for example Fe 1 and Fe 2).
Note (2): All abundances are relative to the Solar values.
     "---"  means that no lines were available to determine the abundance.
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                                      Patricia Vannier [CDS]    31-Mar-2013
