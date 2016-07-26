J/MNRAS/373/13  Lower main-sequence stars fundamental param. (Casagrande+, 2006)
================================================================================
Accurate fundamental parameters for lower main-sequence stars.
    Casagrande L., Portinari L., Flynn C.
   <Mon. Not. R. Astron. Soc., 373, 13-44 (2006)>
   =2006MNRAS.373...13C
================================================================================
ADC_Keywords: Stars, dwarfs ; Stars, early-type ; Photometry ;
              Effective temperatures
Keywords: techniques: photometric - stars: atmospheres -
          stars: fundamental parameters - Hertzsprung-Russell (HR) diagram -
          stars: late-type - infrared: stars

Abstract:
    We derive an empirical effective temperature and bolometric luminosity
    calibration for G and K dwarfs, by applying our own implementation of
    the Infrared Flux Method to multiband photometry. Our study is based
    on 104 stars for which we have excellent BV(RI)_C_ JHK_s_ photometry,
    excellent parallaxes and good metallicities.

    Colours computed from the most recent synthetic libraries (ATLAS9 and
    MARCS) are found to be in good agreement with the empirical colours
    in the optical bands, but some discrepancies still remain in the
    infrared. Synthetic and empirical bolometric corrections also show
    fair agreement.

Description:
    To recover accurate bolometric fluxes and temperatures for the stars,
    we have obtained accurate and homogeneous Johnson-Cousins BV(RI)_C_
    and JHKs photometry for all the 186 stars in our initial sample. For
    most of the stars in our sample with declination north of DE=-25{deg},
    we have made our own photometric observations from April to December
    2004. Observations were done from Finland in full remote mode, using
    the 35-cm telescope piggybacked on the Swedish 60-cm telescope located
    at La Palma in the Canary Islands. A SBIG charge-coupled device was
    used through all the observations. Johnson-Cousins BV(RI)C colours
    were obtained for all stars. Infrared JHKs photometry for the sample
    has been taken from the Two-Micron All-Sky Survey (2MASS) catalogue.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table1.dat    125      104   Observable and physical quantities for our
                              sample stars
--------------------------------------------------------------------------------

See also:
     I/239 : The Hipparcos and Tycho Catalogues (ESA 1997)

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label       Explanations
--------------------------------------------------------------------------------
   1- 10  A10   ---     Name        Name (HD, BD or CD)
  12- 17  I6    ---     HIP         HIP (Cat. I/239) number
  19- 23  F5.2  mas     Plx         Trigonometric parallax
  25- 28  F4.2  mas   e_Plx         rms uncertainty on Plx
  30- 33  I4    K       Teff        Effective temperature (1)
  35- 37  I3    K     e_Teff        rms uncertainty on Teff
  39- 43  F5.3  mas     Diam        ? Angular diameter (1)
  45- 49  F5.3  mas   e_Diam        rms uncertainty on Diam
  51- 56  F6.3  mag     mbol        Apparent bolometric luminosity (2)
  58- 63  F6.3  mag     Vmag        Johnson V magnitude
  65- 69  F5.3  mag     B-V         Johnson B-V colour index
  71- 75  F5.3  mag     V-Rc        Johnson-Cousins V-R colour index
  77- 81  F5.3  mag     R-Ic        Cousins R-I colour index
  83- 86  F4.2  mag     Jmag        2MASS J magnitude
  88- 91  F4.2  mag   e_Jmag        rms uncertainty on Jmag
  93- 96  F4.2  mag     Hmag        2MASS H magnitude
  98-101  F4.2  mag   e_Hmag        rms uncertainty in Hmag
 103-106  F4.2  mag     Ksmag       2MASS Ks magnitude
 108-111  F4.2  mag   e_Ksmag       rms uncertainty on Ks
 113-117  F5.2  [Sun]   [Fe/H]      Metallicity
 119-123  F5.2  [Sun]   [alpha/Fe]  Abundance [alpha/Fe]
     125  A1    ---   r_[alpha/Fe]  Reference for [alpha/Fe] (3)
--------------------------------------------------------------------------------
Note (1): Effective temperatures and angular diameters are those computed
          via IRFM as described in Section 4.
Note (2): Apparent bolometric magnitudes are obtained as described in Section 7,
          where the absolute bolometric magnitude of the Sun is MBol=4.74.
Note (3): Source of metallicities as follows:
      a = Mishenina et al. (2002A&A...396..189M and 2004, Cat. J/A+A/418/551)
      b = Valenti & Fischer (2005)
      c = Pompeia, Barbuy & Grenon (2002ApJ...566..845P + 2003ApJ...592.1173P)
      d = Santos et al. (2005A&A...437.1127S)
      e = Favata, Micela & Sciortino (1997A&A...323..809F)
      f = Gray et al. (2003, Cat. J/AJ/126/2048)
      g = Paulson, Sneden & Cochran (2003AJ....125.3185P)
      h = Santos et al. (2004, Cat. J/A+A/415/1153)
--------------------------------------------------------------------------------

History:
    31-Aug-2007: From electronic version of the journal
    27-Jan-2011: In ReadMe, Plx and Diam columns were inverted

================================================================================
(End)                                      Patricia Vannier [CDS]    08-Jun-2007
