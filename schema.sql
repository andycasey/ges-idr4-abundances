-- ges-idr4-wg15

/*
  THIS IS ONLY RELEVANT FOR THE ABUNDANCE ROUND.
  In a stellar parameter round, the observations would just contain information
  that was not produced by the analysis nodes.
*/
DROP TABLE IF EXISTS observations;
CREATE TABLE observations (
    cname char(21) not null,
    ges_fld char(23) not null,
    object char(28) not null,
    filename char(69) not null,
    ges_type char(8) not null,
    setup char(9) not null,
    wg integer not null,
    ra numeric not null,
    dec numeric not null,
    snr numeric not null,
    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    teff numeric,
    e_teff numeric,
    nn_teff integer,
    enn_teff numeric,
    nne_teff integer,
    sys_err_teff numeric,
    logg numeric,
    e_logg numeric,
    nn_logg integer,
    enn_logg numeric,
    nne_logg integer,
    sys_err_logg numeric,
    feh numeric,
    e_feh numeric,
    nn_feh integer,
    enn_feh numeric,
    nne_feh integer,
    sys_err_feh numeric,
    xi numeric,
    e_xi numeric,
    nn_xi integer,
    enn_xi numeric,
    nne_xi integer,
    mh numeric,
    e_mh numeric,
    nn_mh integer,
    enn_mh numeric,
    nne_mh integer,
    alpha_fe numeric,
    e_alpha_fe numeric,
    nn_alpha_fe integer,
    enn_alpha_fe numeric,
    nne_alpha_fe integer,
    vrad numeric,
    e_vrad numeric,
    vsini numeric,
    e_vsini numeric,
    lim_vsini integer,
    teff_phot numeric,
    e_teff_phot numeric,
    teff_irfm_2mj numeric,
    e_teff_irfm_2mj numeric,
    teff_irfm_2mh numeric,
    e_teff_irfm_2mh numeric,
    teff_irfm_2mks numeric,
    e_teff_irfm_2mks numeric,
    teff_irfm_2mave numeric,
    e_teff_irfm_2mave numeric,
    fbol_irfm_2m numeric);
ALTER TABLE observations ADD COLUMN id BIGSERIAL PRIMARY KEY;


DROP TABLE IF EXISTS flags;
CREATE TABLE flags (
    description char(120) not null);
ALTER TABLE flags ADD COLUMN id BIGSERIAL PRIMARY KEY;


-- Line abundances are provided in ASCII form, primarily from the WG 11 nodes.
DROP TABLE IF EXISTS node_line_abundances;
CREATE TABLE node_line_abundances (
    wg integer not null,
    node char(10) not null,
    cname char(21) not null,
    object char(21) not null,
    abundance_filename char(140) not null,
    spectrum_filename_stub char(140) not null,
    
    measurement_type char(2) not null,
    code char(40) not null,
    
    wavelength numeric not null,
    element char(2) not null,
    ion integer not null,
    loggf numeric default 'NaN',
    chi numeric default 'NaN',

    ew numeric default 'NaN',
    e_ew numeric default 'NaN',
    upper_ew int default 0, 
    
    abundance numeric default 'NaN',
    e_abundance numeric default 'NaN',
    upper_abundance integer default 0,
    flags integer default 0,
    UNIQUE (wg, node, cname, spectrum_filename_stub, element, ion, wavelength)
);
ALTER TABLE node_line_abundances ADD COLUMN id BIGSERIAL PRIMARY KEY;

/*
# The node spectrum abundances come from a given node analysing many lines and
# providing some mean value. These are generally provided in the node result
# FITS files.
*/
DROP TABLE IF EXISTS node_spectrum_abundances;
CREATE TABLE node_spectrum_abundances (
    wg integer not null,
    node char(10) not null,
    cname char(21) not null,
    object char(21) not null,
    filename char(140) not null,
    setup char(9) not null,
    element char(2) not null,
    ion integer not null check (ion > 0),
    abundance numeric default 'NaN',
    e_abundance numeric default 'NaN',
    upper_abundance integer default 0,
    enn numeric default 'NaN',
    nn integer default 0,
    nl integer default 0,
    flags integer default 0,
    UNIQUE (wg, node, filename, element, ion)
);
ALTER TABLE node_spectrum_abundances ADD COLUMN id BIGSERIAL PRIMARY KEY;



-- The WG-level abundances are unique on a per (WG, filename, element, ion) level.
DROP TABLE IF EXISTS wg_abundances;
CREATE TABLE wg_abundances (
    wg integer not null,
    cname char(21) not null,
    filename char(140) not null,
    element char(2) not null,
    ion integer not null check (ion > 0),
    abundance numeric default 'NaN',
    e_abundance numeric default 'NaN',
    enn numeric default 'NaN',
    nn integer default 0,
    nl integer default 0,
    upper_abundance integer default 0,
    flags integer default 0,
    UNIQUE (wg, filename, element, ion)
);
ALTER TABLE wg_abundances ADD COLUMN id BIGSERIAL PRIMARY KEY;



-- The Survey-level abundances are unique on a per (CNAME, element, ion) level.
DROP TABLE IF EXISTS survey_abundances;
CREATE TABLE survey_abundances (
    id_provenance integer not null,
    wg_provenance integer not null,
    cname char(21) not null,
    element char(2) not null,
    ion integer not null check (ion > 0),
    abundance numeric default 'NaN',
    e_abundance numeric default 'NaN',
    upper_abundance numeric default 0,
    enn numeric default 'NaN',
    nn integer default 0,
    nl integer default 0,
    flags integer default 0,
    UNIQUE (cname, element, ion)
);
ALTER TABLE survey_abundances ADD COLUMN id BIGSERIAL PRIMARY KEY;

