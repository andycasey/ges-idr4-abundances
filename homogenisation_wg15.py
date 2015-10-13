

""" Script to remove spurious WG12 recommended measurements. """

import release

ges = release.DataRelease("arc", password="password", host="/tmp")


# The following species have very different results with respect to stars in
# common and WG11.
ges.execute("""UPDATE wg_homogenised_abundances
    SET wg15_abundance = 'NaN', flag = '1'
    WHERE wg = '12' AND (
        (element = 'Mn' AND ion = '1') OR
        (element = 'Co' AND ion = '1') OR
        (element = 'Si' AND ion = '2') OR
        (trim(element) = 'V' AND ion = '1')
    )""")

# La 2 abundances are bad for dwarfs (logg > 3)
ges.execute("""UPDATE wg_homogenised_abundances
    SET wg15_abundance = 'NaN', flag = '1'
    WHERE wg = '12' AND element = 'La' and ion = '2'
        AND cname in (
            SELECT DISTINCT ON (cname) cname
            FROM observations WHERE logg > 3)
    """)
ges.commit()